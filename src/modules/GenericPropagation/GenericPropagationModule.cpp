/**
 * @file
 * @brief Implementation of generic charge propagation module
 * @remarks Based on code from Paul Schuetze
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "GenericPropagationModule.hpp"

#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <random>
#include <sstream>
#include <string>
#include <utility>

#include <Eigen/Core>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TPaveText.h>
#include <TPolyLine3D.h>
#include <TPolyMarker3D.h>
#include <TStyle.h>

#include "core/config/Configuration.hpp"
#include "core/messenger/Messenger.hpp"
#include "core/utils/file.h"
#include "core/utils/log.h"
#include "core/utils/unit.h"
#include "tools/ROOT.h"
#include "tools/runge_kutta.h"

#include "objects/DepositedCharge.hpp"
#include "objects/PropagatedCharge.hpp"

using namespace allpix;

/**
 * Besides binding the message and setting defaults for the configuration, the module copies some configuration variables to
 * local copies to speed up computation.
 */
GenericPropagationModule::GenericPropagationModule(Configuration config,
                                                   Messenger* messenger,
                                                   std::shared_ptr<Detector> detector)
    : Module(std::move(config), detector), messenger_(messenger), detector_(std::move(detector)) {
    // Enable parallelization of this module if multithreading is enabled
    enable_parallelization();

    // Save detector model
    model_ = detector_->getModel();

    // Require deposits message for single detector
    messenger_->bindSingle(this, &GenericPropagationModule::deposits_message_, MsgFlags::REQUIRED);

    // Seed the random generator with the module seed
    random_generator_.seed(getRandomSeed());

    // Set default value for config variables
    config_.setDefault<double>("spatial_precision", Units::get(0.25, "nm"));
    config_.setDefault<double>("timestep_start", Units::get(0.01, "ns"));
    config_.setDefault<double>("timestep_min", Units::get(0.001, "ns"));
    config_.setDefault<double>("timestep_max", Units::get(0.5, "ns"));
    config_.setDefault<double>("integration_time", Units::get(25, "ns"));
    config_.setDefault<unsigned int>("charge_per_step", 10);
    config_.setDefault<double>("temperature", 293.15);

    config_.setDefault<bool>("output_plots", false);
    config_.setDefault<bool>("output_animations", false);
    config_.setDefault<bool>("output_animations_color_markers", false);
    config_.setDefault<double>("output_plots_step", config_.get<double>("timestep_max"));
    config_.setDefault<bool>("output_plots_use_pixel_units", false);
    config_.setDefault<bool>("output_plots_align_pixels", false);
    config_.setDefault<double>("output_plots_theta", 0.0f);
    config_.setDefault<double>("output_plots_phi", 0.0f);

    // Set defaults for charge carrier propagation:
    config_.setDefault<bool>("propagate_electrons", true);
    config_.setDefault<bool>("propagate_holes", false);
    if(!config_.get<bool>("propagate_electrons") && !config_.get<bool>("propagate_holes")) {
        throw InvalidValueError(
            config_,
            "propagate_electrons",
            "No charge carriers selected for propagation, enable 'propagate_electrons' or 'propagate_holes'.");
    }
    // Copy some variables from configuration to avoid lookups:
    temperature_ = config_.get<double>("temperature");
    timestep_min_ = config_.get<double>("timestep_min");
    timestep_max_ = config_.get<double>("timestep_max");
    timestep_start_ = config_.get<double>("timestep_start");
    integration_time_ = config_.get<double>("integration_time");
    target_spatial_precision_ = config_.get<double>("spatial_precision");
    output_plots_ = config_.get<bool>("output_plots");
    output_plots_step_ = config_.get<double>("output_plots_step");

    // Parameterization variables from https://doi.org/10.1016/0038-1101(77)90054-5 (section 5.2)
    electron_Vm_ = Units::get(1.53e9 * std::pow(temperature_, -0.87), "cm/s");
    electron_Ec_ = Units::get(1.01 * std::pow(temperature_, 1.55), "V/cm");
    electron_Beta_ = 2.57e-2 * std::pow(temperature_, 0.66);

    hole_Vm_ = Units::get(1.62e8 * std::pow(temperature_, -0.52), "cm/s");
    hole_Ec_ = Units::get(1.24 * std::pow(temperature_, 1.68), "V/cm");
    hole_Beta_ = 0.46 * std::pow(temperature_, 0.17);

    boltzmann_kT_ = Units::get(8.6173e-5, "eV/K") * temperature_;
    elementary_charge_ = 1.6021766208e-19; // C = As
}

void GenericPropagationModule::create_output_plots(unsigned int event_num) {
    LOG(TRACE) << "Writing output plots";

    // Convert to pixel units if necessary
    if(config_.get<bool>("output_plots_use_pixel_units")) {
        for(auto& deposit_points : output_plot_points_) {
            for(auto& point : deposit_points.second) {
                point.SetX((point.x() / model_->getPixelSize().x()) + 1);
                point.SetY((point.y() / model_->getPixelSize().y()) + 1);
            }
        }
    }

    // Calculate the axis limits
    double minX = FLT_MAX, maxX = FLT_MIN;
    double minY = FLT_MAX, maxY = FLT_MIN;
    unsigned long tot_point_cnt = 0;
    double start_time = std::numeric_limits<double>::max();
    unsigned int total_charge = 0;
    unsigned int max_charge = 0;
    for(auto& deposit_points : output_plot_points_) {
        for(auto& point : deposit_points.second) {
            minX = std::min(minX, point.x());
            maxX = std::max(maxX, point.x());

            minY = std::min(minY, point.y());
            maxY = std::max(maxY, point.y());
        }
        start_time = std::min(start_time, deposit_points.first.getEventTime());
        total_charge += deposit_points.first.getCharge();
        max_charge = std::max(max_charge, deposit_points.first.getCharge());

        tot_point_cnt += deposit_points.second.size();
    }

    // Compute frame axis sizes if equal scaling is requested
    if(config_.get<bool>("output_plots_use_equal_scaling", true)) {
        double centerX = (minX + maxX) / 2.0;
        double centerY = (minY + maxY) / 2.0;
        if(config_.get<bool>("output_plots_use_pixel_units")) {
            minX = centerX - model_->getSensorSize().z() / model_->getPixelSize().x() / 2.0;
            maxX = centerX + model_->getSensorSize().z() / model_->getPixelSize().x() / 2.0;

            minY = centerY - model_->getSensorSize().z() / model_->getPixelSize().y() / 2.0;
            maxY = centerY + model_->getSensorSize().z() / model_->getPixelSize().y() / 2.0;
        } else {
            minX = centerX - model_->getSensorSize().z() / 2.0;
            maxX = centerX + model_->getSensorSize().z() / 2.0;

            minY = centerY - model_->getSensorSize().z() / 2.0;
            maxY = centerY + model_->getSensorSize().z() / 2.0;
        }
    }

    // Align on pixels if requested
    if(config_.get<bool>("output_plots_align_pixels")) {
        if(config_.get<bool>("output_plots_use_pixel_units")) {
            minX = std::floor(minX - 0.5) + 0.5;
            minY = std::floor(minY + 0.5) - 0.5;
            maxX = std::ceil(maxX - 0.5) + 0.5;
            maxY = std::ceil(maxY + 0.5) - 0.5;
        } else {
            double div;
            div = minX / model_->getPixelSize().x();
            minX = (std::floor(div - 0.5) + 0.5) * model_->getPixelSize().x();
            div = minY / model_->getPixelSize().y();
            minY = (std::floor(div - 0.5) + 0.5) * model_->getPixelSize().y();
            div = maxX / model_->getPixelSize().x();
            maxX = (std::ceil(div + 0.5) - 0.5) * model_->getPixelSize().x();
            div = maxY / model_->getPixelSize().y();
            maxY = (std::ceil(div + 0.5) - 0.5) * model_->getPixelSize().y();
        }
    }

    // Create 2D histogram for deposition positions
    // Find axis range
    double minXDep = FLT_MAX; 
    double maxXDep = FLT_MIN;
    double minYDep = FLT_MAX;
    double maxYDep = FLT_MIN;
    for(auto& position : output_plot_points_initial_) {
        minXDep = std::min(minXDep, position.x());
        maxXDep = std::max(maxXDep, position.x());

        minYDep = std::min(minYDep, position.y());
        maxYDep = std::max(maxYDep, position.y());
    }

    auto canvas_deposition = std::make_unique<TCanvas>(("deposition_plot_" + std::to_string(event_num)).c_str(),
                                            ("Deposition position projection onto pixel layer for event " + std::to_string(event_num)).c_str(),
                                            1280,
                                            1024);
    canvas_deposition->cd();
    auto* histogram_deposition = new TH2F(("deposition_plot_" + std::to_string(event_num)).c_str(),
                                            ("Deposition position projection onto pixel layer for event " + std::to_string(event_num)).c_str(),
                                            200, minXDep, maxXDep,
                                            200, minYDep, maxYDep);

    for(auto& position : output_plot_points_initial_) {
        histogram_deposition->Fill(position.x(), position.y());
    }

    histogram_deposition->Draw("COLZ");
    histogram_deposition->Write();
    getROOTDirectory()->WriteTObject(canvas_deposition.get());

    // Create 2D histogram for charge projection onto pixel layer
    auto canvas_pixel = std::make_unique<TCanvas>(("pixel_plot_" + std::to_string(event_num)).c_str(),
                                            ("Charge projection onto pixel layer for event " + std::to_string(event_num)).c_str(),
                                            1280,
                                            1024);
    canvas_pixel->cd();
    auto* histogram_pixel = new TH2F(("pixel_plot_" + std::to_string(event_num)).c_str(),
                                            ("Charge projection onto pixel layer for event " + std::to_string(event_num)).c_str(),
                                            200, minX, maxX,
                                            200, minY, maxY);

    for(auto& position : output_plot_points_final_) {
        histogram_pixel->Fill(position.x(), position.y());
    }

    histogram_pixel->Draw("COLZ");
    histogram_pixel->Write();
    getROOTDirectory()->WriteTObject(canvas_pixel.get());

    // Use a histogram to create the underlying frame
    auto* histogram_frame = new TH3F(("frame_" + getUniqueName() + "_" + std::to_string(event_num)).c_str(),
                                     "",
                                     10,
                                     minX,
                                     maxX,
                                     10,
                                     minY,
                                     maxY,
                                     10,
                                     model_->getSensorCenter().z() - model_->getSensorSize().z() / 2.0,
                                     model_->getSensorCenter().z() + model_->getSensorSize().z() / 2.0);
    histogram_frame->SetDirectory(getROOTDirectory());

    // Create the canvas for the line plot and set orientation
    auto canvas = std::make_unique<TCanvas>(("line_plot_" + std::to_string(event_num)).c_str(),
                                            ("Propagation of charge for event " + std::to_string(event_num)).c_str(),
                                            1280,
                                            1024);
    canvas->cd();
    canvas->SetTheta(config_.get<float>("output_plots_theta") * 180.0f / ROOT::Math::Pi());
    canvas->SetPhi(config_.get<float>("output_plots_phi") * 180.0f / ROOT::Math::Pi());

    // Draw the frame on the canvas
    histogram_frame->GetXaxis()->SetTitle(
        (std::string("x ") + (config_.get<bool>("output_plots_use_pixel_units") ? "(pixels)" : "(mm)")).c_str());
    histogram_frame->GetYaxis()->SetTitle(
        (std::string("y ") + (config_.get<bool>("output_plots_use_pixel_units") ? "(pixels)" : "(mm)")).c_str());
    histogram_frame->GetZaxis()->SetTitle("z (mm)");
    histogram_frame->Draw();

    // Loop over all point sets created during propagation
    // The vector of unique_pointers is required in order not to delete the objects before the canvas is drawn.
    std::vector<std::unique_ptr<TPolyLine3D>> lines;
    short current_color = 1;
    for(auto& deposit_points : output_plot_points_) {
        auto line = std::make_unique<TPolyLine3D>();
        for(auto& point : deposit_points.second) {
            line->SetNextPoint(point.x(), point.y(), point.z());
        }
        // Plot all lines with at least three points with different color
        if(line->GetN() >= 3) {
            EColor plot_color = (deposit_points.first.getType() == CarrierType::ELECTRON ? EColor::kAzure : EColor::kOrange);
            current_color = static_cast<short int>(plot_color - 9 + (static_cast<int>(current_color) + 1) % 19);
            line->SetLineColor(current_color);
            line->Draw("same");
        }
        lines.push_back(std::move(line));
    }

    // Draw and write canvas to module output file, then clear the stored lines
    canvas->Draw();
    getROOTDirectory()->WriteTObject(canvas.get());
    lines.clear();

    // Create canvas for GIF animition of process
    canvas = std::make_unique<TCanvas>(("animation_" + std::to_string(event_num)).c_str(),
                                       ("Propagation of charge for event " + std::to_string(event_num)).c_str(),
                                       1280,
                                       1024);
    canvas->cd();

    // Change axis labels if close to zero or PI as they behave different here
    if(std::fabs(config_.get<double>("output_plots_theta") / (ROOT::Math::Pi() / 2.0) -
                 std::round(config_.get<double>("output_plots_theta") / (ROOT::Math::Pi() / 2.0))) < 1e-6 ||
       std::fabs(config_.get<double>("output_plots_phi") / (ROOT::Math::Pi() / 2.0) -
                 std::round(config_.get<double>("output_plots_phi") / (ROOT::Math::Pi() / 2.0))) < 1e-6) {
        histogram_frame->GetXaxis()->SetLabelOffset(-0.1f);
        histogram_frame->GetYaxis()->SetLabelOffset(-0.075f);
    } else {
        histogram_frame->GetXaxis()->SetTitleOffset(2.0f);
        histogram_frame->GetYaxis()->SetTitleOffset(2.0f);
    }

    // Draw frame on canvas
    histogram_frame->Draw();

    if(config_.get<bool>("output_animations")) {
        // Create the contour histogram
        std::vector<std::string> file_name_contour;
        std::vector<TH2F*> histogram_contour;
        file_name_contour.push_back(createOutputFile("contourX" + std::to_string(event_num) + ".gif"));
        histogram_contour.push_back(new TH2F(("contourX_" + getUniqueName() + "_" + std::to_string(event_num)).c_str(),
                                             "",
                                             100,
                                             minY,
                                             maxY,
                                             100,
                                             model_->getSensorCenter().z() - model_->getSensorSize().z() / 2.0,
                                             model_->getSensorCenter().z() + model_->getSensorSize().z() / 2.0));
        histogram_contour.back()->SetDirectory(getROOTDirectory());
        file_name_contour.push_back(createOutputFile("contourY" + std::to_string(event_num) + ".gif"));
        histogram_contour.push_back(new TH2F(("contourY_" + getUniqueName() + "_" + std::to_string(event_num)).c_str(),
                                             "",
                                             100,
                                             minX,
                                             maxX,
                                             100,
                                             model_->getSensorCenter().z() - model_->getSensorSize().z() / 2.0,
                                             model_->getSensorCenter().z() + model_->getSensorSize().z() / 2.0));
        histogram_contour.back()->SetDirectory(getROOTDirectory());
        file_name_contour.push_back(createOutputFile("contourZ" + std::to_string(event_num) + ".gif"));
        histogram_contour.push_back(new TH2F(("contourZ_" + getUniqueName() + "_" + std::to_string(event_num)).c_str(),
                                             "",
                                             100,
                                             minX,
                                             maxX,
                                             100,
                                             minY,
                                             maxY));
        histogram_contour.back()->SetDirectory(getROOTDirectory());

        // Create file and disable statistics for histogram
        std::string file_name_anim = createOutputFile("animation" + std::to_string(event_num) + ".gif");
        for(size_t i = 0; i < 3; ++i) {
            histogram_contour[i]->SetStats(false);
        }

        // Remove temporary created files
        remove(file_name_anim.c_str());
        for(size_t i = 0; i < 3; ++i) {
            remove(file_name_contour[i].c_str());
        }

        // Create color table
        TColor* colors[80];
        for(int i = 20; i < 100; ++i) {
            auto color_idx = TColor::GetFreeColorIndex();
            colors[i - 20] = new TColor(color_idx,
                                        static_cast<float>(i) / 100.0f - 0.2f,
                                        static_cast<float>(i) / 100.0f - 0.2f,
                                        static_cast<float>(i) / 100.0f - 0.2f);
        }

        // Create animation of moving charges
        auto animation_time = static_cast<unsigned int>(
            std::round((Units::convert(config_.get<long double>("output_plots_step"), "ms") / 10.0) *
                       config_.get<long double>("output_animations_time_scaling", 1e9)));
        unsigned long plot_idx = 0;
        unsigned int point_cnt = 0;
        LOG_PROGRESS(INFO, getUniqueName() + "_OUTPUT_PLOTS") << "Written 0 of " << tot_point_cnt << " points for animation";
        while(point_cnt < tot_point_cnt) {
            std::vector<std::unique_ptr<TPolyMarker3D>> markers;
            unsigned long min_idx_diff = std::numeric_limits<unsigned long>::max();

            // Reset the canvas
            canvas->Clear();
            canvas->SetTheta(config_.get<float>("output_plots_theta") * 180.0f / ROOT::Math::Pi());
            canvas->SetPhi(config_.get<float>("output_plots_phi") * 180.0f / ROOT::Math::Pi());
            canvas->Draw();

            // Reset the histogram frame
            histogram_frame->SetTitle("Charge propagation in sensor");
            histogram_frame->GetXaxis()->SetTitle(
                (std::string("x ") + (config_.get<bool>("output_plots_use_pixel_units") ? "(pixels)" : "(mm)")).c_str());
            histogram_frame->GetYaxis()->SetTitle(
                (std::string("y ") + (config_.get<bool>("output_plots_use_pixel_units") ? "(pixels)" : "(mm)")).c_str());
            histogram_frame->GetZaxis()->SetTitle("z (mm)");
            histogram_frame->Draw();

            auto text = std::make_unique<TPaveText>(-0.75, -0.75, -0.60, -0.65);
            auto time_ns = Units::convert(plot_idx * config_.get<long double>("output_plots_step"), "ns");
            std::stringstream sstr;
            sstr << std::fixed << std::setprecision(2) << time_ns << "ns";
            auto time_str = std::string(8 - sstr.str().size(), ' ');
            time_str += sstr.str();
            text->AddText(time_str.c_str());
            text->Draw();

            // Plot all the required points
            for(auto& deposit_points : output_plot_points_) {
                auto points = deposit_points.second;

                auto diff = static_cast<unsigned long>(std::round((deposit_points.first.getEventTime() - start_time) /
                                                                  config_.get<long double>("output_plots_step")));
                if(static_cast<long>(plot_idx) - static_cast<long>(diff) < 0) {
                    min_idx_diff = std::min(min_idx_diff, diff - plot_idx);
                    continue;
                }
                auto idx = plot_idx - diff;
                if(idx >= points.size()) {
                    continue;
                }
                min_idx_diff = 0;

                auto marker = std::make_unique<TPolyMarker3D>();
                marker->SetMarkerStyle(kFullCircle);
                marker->SetMarkerSize(static_cast<float>(deposit_points.first.getCharge() *
                                                         config_.get<unsigned int>("output_animations_marker_size", 1)) /
                                      static_cast<float>(max_charge));
                auto initial_z_perc = static_cast<int>(
                    ((points[0].z() + model_->getSensorSize().z() / 2.0) / model_->getSensorSize().z()) * 80);
                initial_z_perc = std::max(std::min(79, initial_z_perc), 0);
                if(config_.get<bool>("output_animations_color_markers")) {
                    marker->SetMarkerColor(static_cast<Color_t>(colors[initial_z_perc]->GetNumber()));
                }
                marker->SetNextPoint(points[idx].x(), points[idx].y(), points[idx].z());
                marker->Draw();
                markers.push_back(std::move(marker));

                histogram_contour[0]->Fill(points[idx].y(), points[idx].z(), deposit_points.first.getCharge());
                histogram_contour[1]->Fill(points[idx].x(), points[idx].z(), deposit_points.first.getCharge());
                histogram_contour[2]->Fill(points[idx].x(), points[idx].y(), deposit_points.first.getCharge());
                ++point_cnt;
            }

            // Create a step in the animation
            if(min_idx_diff != 0) {
                canvas->Print((file_name_anim + "+100").c_str());
                plot_idx += min_idx_diff;
            } else {
                // print animation
                if(point_cnt < tot_point_cnt - 1) {
                    canvas->Print((file_name_anim + "+" + std::to_string(animation_time)).c_str());
                } else {
                    canvas->Print((file_name_anim + "++100").c_str());
                }

                // Draw and print contour histograms
                for(size_t i = 0; i < 3; ++i) {
                    canvas->Clear();
                    canvas->SetTitle((std::string("Contour of charge propagation projected on the ") +
                                      static_cast<char>('X' + i) + "-axis")
                                         .c_str());
                    switch(i) {
                    case 0 /* x */:
                        histogram_contour[i]->GetXaxis()->SetTitle(
                            (std::string("y ") + (config_.get<bool>("output_plots_use_pixel_units") ? "(pixels)" : "(mm)"))
                                .c_str());
                        histogram_contour[i]->GetYaxis()->SetTitle("z (mm)");
                        break;
                    case 1 /* y */:
                        histogram_contour[i]->GetXaxis()->SetTitle(
                            (std::string("x ") + (config_.get<bool>("output_plots_use_pixel_units") ? "(pixels)" : "(mm)"))
                                .c_str());
                        histogram_contour[i]->GetYaxis()->SetTitle("z (mm)");
                        break;
                    case 2 /* z */:
                        histogram_contour[i]->GetXaxis()->SetTitle(
                            (std::string("x ") + (config_.get<bool>("output_plots_use_pixel_units") ? "(pixels)" : "(mm)"))
                                .c_str());
                        histogram_contour[i]->GetYaxis()->SetTitle(
                            (std::string("y ") + (config_.get<bool>("output_plots_use_pixel_units") ? "(pixels)" : "(mm)"))
                                .c_str());
                        break;
                    default:;
                    }
                    histogram_contour[i]->SetMinimum(1);
                    histogram_contour[i]->SetMaximum(total_charge /
                                                     config_.get<double>("output_plots_contour_max_scaling", 10));
                    histogram_contour[i]->Draw("CONTZ 0");
                    if(point_cnt < tot_point_cnt - 1) {
                        canvas->Print((file_name_contour[i] + "+" + std::to_string(animation_time)).c_str());
                    } else {
                        canvas->Print((file_name_contour[i] + "++100").c_str());
                    }
                    histogram_contour[i]->Reset();
                }
                ++plot_idx;
            }
            markers.clear();

            LOG_PROGRESS(INFO, getUniqueName() + "_OUTPUT_PLOTS")
                << "Written " << point_cnt << " of " << tot_point_cnt << " points for animation";
        }
    }
    output_plot_points_.clear();
    output_plot_points_final_.clear();
    output_plot_points_initial_.clear();
}

void GenericPropagationModule::init() {
    auto detector = getDetector();

    // Check for electric field and output warning for slow propagation if not defined
    if(!detector->hasElectricField()) {
        LOG(WARNING) << "This detector does not have an electric field.";
    }

    // For linear fields we can in addition check if the correct carriers are propagated
    if(detector->getElectricFieldType() == ElectricFieldType::LINEAR) {
        auto model = detector_->getModel();
        auto probe_point = ROOT::Math::XYZPoint(model->getSensorCenter().x(),
                                                model->getSensorCenter().y(),
                                                model->getSensorCenter().z() + model->getSensorSize().z() / 2.01);

        // Get the field close to the implants and check its sign:
        auto efield = detector->getElectricField(probe_point);
        auto direction = std::signbit(efield.z());
        // Compare with propagated carrier type:
        if(direction && !config_.get<bool>("propagate_electrons")) {
            LOG(WARNING) << "Electric field indicates electron collection at implants, but electrons are not propagated!";
        }
        if(!direction && !config_.get<bool>("propagate_holes")) {
            LOG(WARNING) << "Electric field indicates hole collection at implants, but holes are not propagated!";
        }
    }
    if(output_plots_) {
        auto time_bins =
            static_cast<int>(Units::convert(integration_time_ / config_.get<long double>("output_plots_step"), "ns"));
        drift_time_histo = new TH1D(
            "drift_time_histo", "Drift time;t[ns];charge carriers", time_bins, 0., static_cast<int>(integration_time_));
    }
}

void GenericPropagationModule::run(unsigned int event_num) {

    // Create vector of propagated charges to output
    std::vector<PropagatedCharge> propagated_charges;

    // Loop over all deposits for propagation
    LOG(TRACE) << "Propagating charges in sensor";
    unsigned int propagated_charges_count = 0;
    unsigned int step_count = 0;
    long double total_time = 0;

    // unsigned int charge_per_step = config_.get<unsigned int>("charge_per_step");
    // Get length of track
    // TODO: Fix - ugly!
    /*
    double track_length_ = 0;
    double total_charge_ = 0;
    auto& dp = deposits_message_->getData();
    auto initial_pos = dp.at(0).getLocalPosition();

    // std::advance(dp, 1);
    for(auto& deposit : dp) {
        auto position = deposit.getLocalPosition();
        double step_length = std::sqrt((position - initial_pos).Mag2());
        // LOG(INFO) << step_length << " " << deposit.getCharge();
        if(step_length < Units::get(10., "um")) {
            track_length_ += step_length;
        }
        initial_pos = position;

        total_charge_ += deposit.getCharge();
    }


    // Get the charge_density of the currently processed track
    charge_density_ = 0.5 * total_charge_ * elementary_charge_ / track_length_;
    //LOG(INFO) << "Track length: " << track_length_;
    // LOG(INFO) << "Energy: " << 0.5 * total_charge_ * 3.64 / 1000;
    // LOG(INFO) << "Ratio: " << 1./ (1 + 2 * boltzmann_kT_ / charge_density_  * ROOT::Math::Pi() * 11.68 * 8.854187817e-12 * 10000 );
    // LOG(INFO) << "Ratio: " << charge_density_ / (charge_density_ + 2*boltzmann_kT_*ROOT::Math::Pi()*11.68 * 8.854187817e-12 * 1000);
    */ 

    // auto position_first = deposits_message_->getData().at(0).getLocalPosition();
    // auto position_second = deposits_message_->getData().at(0).getLocalPosition();

    // for(auto& deposit : deposits_message_->getData()) {
    bool cylinder;
    for(auto it = deposits_message_->getData().begin(); it != deposits_message_->getData().end(); it++) {
        auto& deposit = *it;
        if((deposit.getType() == CarrierType::ELECTRON && !config_.get<bool>("propagate_electrons")) ||
           (deposit.getType() == CarrierType::HOLE && !config_.get<bool>("propagate_holes"))) {
            LOG(DEBUG) << "Skipping charge carriers (" << deposit.getType() << ") on "
                       << display_vector(deposit.getLocalPosition(), {"mm", "um"});
            continue;
        }

        // Loop over all charges in the deposit
        unsigned int charges_remaining = deposit.getCharge();
        // LOG(INFO) << "Step: " << step_count << ", Charge: " << charges_remaining;

        LOG(DEBUG) << "Set of charge carriers (" << deposit.getType() << ") on "
                   << display_vector(deposit.getLocalPosition(), {"mm", "um"});

        // Get position and propagate through sensor
        // auto position_first = deposit.getLocalPosition();
        // std::advance(deposit);
        auto position_first = deposit.getLocalPosition();
        ROOT::Math::XYZVector step_direction;
        double step_length = Units::get(20., "um");

        if(std::next(it, 2) != deposits_message_->getData().end()) {
            auto& second_deposit = *std::next(it, 2);
            auto position_second = second_deposit.getLocalPosition();
            step_direction = position_second - position_first;
            // LOG(INFO) << position_first;
            // LOG(INFO) << position_second;

            // unsigned int charge_second = deposit.getCharge();
            step_length = std::sqrt((position_second - position_first).Mag2());
        }
        // LOG(INFO) << "Step length: " << step_length << ", Limit: " << Units::get(10., "um");
        // LOG(INFO) << Units::get(step_length, "um");

        if(step_length < Units::get(1., "um")) {
            v_vec = {step_direction.x(), step_direction.y(), step_direction.z()};
            charge_density_ = charges_remaining * elementary_charge_ / step_length;
            // LOG(INFO) << charge_density_;
            cylinder = true;
        } else {
            // LOG(INFO) << "Sphere";
            step_direction = ROOT::Math::XYZVector(0, 0, 0);
            v_vec = {0, 0, 0};
            charge_density_ = 3. / 4 * elementary_charge_ * charges_remaining;
            cylinder = false;
        }
        // LOG(INFO) << "Cylinder? " << cylinder;

        auto charge_per_step = config_.get<unsigned int>("charge_per_step");
        double deposit_total_steps = static_cast<int>(charges_remaining / charge_per_step + 0.5);
        while(charges_remaining > 0) {
            // Define number of charges to be propagated and remove charges of this step from the total
            if(charge_per_step > charges_remaining) {
                charge_per_step = charges_remaining;
            }
            charges_remaining -= charge_per_step;

            // Add point of deposition to the output plots if requested
            if(output_plots_) {
                auto global_position = detector_->getGlobalPosition(position_first);
                output_plot_points_.emplace_back(
                    PropagatedCharge(position_first, global_position, deposit.getType(), charge_per_step, deposit.getEventTime()),
                    std::vector<ROOT::Math::XYZPoint>());
        
                // For 2D histogram of initial deposition points
                output_plot_points_initial_.push_back(position_first);
            }

            // Propagate a single charge deposit
            auto prop_pair = propagate(position_first, deposit.getType(), cylinder);
            auto position = prop_pair.first;

            // Advance one step
            position_first += 1./deposit_total_steps * step_direction;

            LOG(DEBUG) << " Propagated " << charge_per_step << " to " << display_vector(position, {"mm", "um"}) << " in "
                       << Units::display(prop_pair.second, "ns") << " time";

            // Create a new propagated charge and add it to the list
            auto global_position = detector_->getGlobalPosition(position);
            PropagatedCharge propagated_charge(position,
                                               global_position,
                                               deposit.getType(),
                                               charge_per_step,
                                               deposit.getEventTime() + prop_pair.second,
                                               &deposit);

            propagated_charges.push_back(std::move(propagated_charge));
            if(output_plots_) {
                output_plot_points_final_.push_back(position);
            }

            // Update statistical information
            ++step_count;
            propagated_charges_count += charge_per_step;
            total_time += charge_per_step * prop_pair.second;

            // Fill plot for drift time
            if(output_plots_) {
                drift_time_histo->SetBinContent(
                    drift_time_histo->FindBin(prop_pair.second),
                    drift_time_histo->GetBinContent(drift_time_histo->FindBin(prop_pair.second)) + charge_per_step);
            }
        }

        // Get next position of track segment
        // position_first = position_second;
    }

    // Output plots if required
    if(output_plots_) {
        create_output_plots(event_num);
    }

    // Write summary and update statistics
    long double average_time = total_time / std::max(1u, propagated_charges_count);
    LOG(INFO) << "Propagated " << propagated_charges_count << " charges in " << step_count << " steps in average time of "
              << Units::display(average_time, "ns");
    total_propagated_charges_ += propagated_charges_count;
    total_steps_ += step_count;
    total_time_ += total_time;

    // Create a new message with propagated charges
    auto propagated_charge_message = std::make_shared<PropagatedChargeMessage>(std::move(propagated_charges), detector_);

    // Dispatch the message with propagated charges
    messenger_->dispatchMessage(this, propagated_charge_message);
}

/**
 * Propagation is simulated using a parameterization for the electron mobility. This is used to calculate the electron
 * velocity at every point with help of the electric field map of the detector. An Runge-Kutta integration is applied in
 * multiple steps, adding a random diffusion to the propagating charge every step.
 */
std::pair<ROOT::Math::XYZPoint, double> GenericPropagationModule::propagate(const ROOT::Math::XYZPoint& pos,
                                                                            const CarrierType& type,
                                                                            bool cylinder) {
    // Create a runge kutta solver using the electric field as step function
    Eigen::Vector3d position(pos.x(), pos.y(), pos.z());

    // Define a lambda function to compute the carrier mobility
    // NOTE This function is typically the most frequently executed part of the framework and therefore the bottleneck
    auto carrier_mobility = [&](double efield_mag) {
        // Compute carrier mobility from constants and electric field magnitude
        double numerator, denominator;
        if(type == CarrierType::ELECTRON) {
            numerator = electron_Vm_ / electron_Ec_;
            denominator = std::pow(1. + std::pow(efield_mag / electron_Ec_, electron_Beta_), 1.0 / electron_Beta_);
        } else {
            numerator = hole_Vm_ / hole_Ec_;
            denominator = std::pow(1. + std::pow(efield_mag / hole_Ec_, hole_Beta_), 1.0 / hole_Beta_);
        }
        return numerator / denominator;
    };

    // Define a function to compute the diffusion
    auto carrier_diffusion = [&](double efield_mag, double timestep) -> Eigen::Vector3d {
        double diffusion_constant = boltzmann_kT_ * carrier_mobility(efield_mag);
        double diffusion_std_dev = std::sqrt(2. * diffusion_constant * timestep);
        // double diffusion_coeff = 2. * diffusion_constant * timestep;

        // double diffusion_std_dev = std::sqrt(diffusion_coeff + repulsion_coeff);
        // LOG(INFO) << "Diffusion " << diffusion_std_dev;
        // LOG(INFO) << "Repulsion " << repulsion_std_dev;
        // LOG(INFO) << "";

        // Get orthogonal vectors to the normal vector v
        // std::vector<double> v_vec = std::vector<double>(vector_v.X(), vector_v.Y(), vector_v.Z());

        // Compute the independent diffusion in three dimensions
        std::normal_distribution<double> gauss_distribution(0, diffusion_std_dev);

        Eigen::Vector3d diffusion;
        for(int i = 0; i < 3; ++i) {
            diffusion[i] = gauss_distribution(random_generator_);
        }

        // LOG(INFO) << "Repulsion: " << repulsion.norm();
        // LOG(INFO) << "Diffusion: " << diffusion.norm();
        return diffusion; //  + repulsion;
    };

    // Define a lambda function to compute the electron velocity
    auto carrier_velocity = [&](double, Eigen::Vector3d cur_pos) -> Eigen::Vector3d {
        auto raw_field = detector_->getElectricField(static_cast<ROOT::Math::XYZPoint>(cur_pos));
        // Compute the drift velocity
        Eigen::Vector3d efield(raw_field.x(), raw_field.y(), raw_field.z());

        return static_cast<int>(type) * carrier_mobility(efield.norm()) * efield;
    };

    // Create the runge kutta solver with an RKF5 tableau
    auto runge_kutta = make_runge_kutta(tableau::RK5, carrier_velocity, timestep_start_, position);

    // Continue propagation until the deposit is outside the sensor
    Eigen::Vector3d last_position = position;
    double last_time = 0;
    size_t next_idx = 0;

    // Get initial mobility
    double carrier_mobility_init = carrier_mobility(std::sqrt((detector_->getElectricField(static_cast<ROOT::Math::XYZPoint>(position)).Mag2())));

    while(detector_->isWithinSensor(static_cast<ROOT::Math::XYZPoint>(position)) &&
          runge_kutta.getTime() < integration_time_) {
        // Update output plots if necessary (depending on the plot step)
        if(output_plots_) {
            auto time_idx = static_cast<size_t>(runge_kutta.getTime() / output_plots_step_);
            while(next_idx <= time_idx) {
                output_plot_points_.back().second.push_back(static_cast<ROOT::Math::XYZPoint>(position));
                next_idx = output_plot_points_.back().second.size();
            }
        }

        // Save previous position and time
        last_position = position;
        last_time = runge_kutta.getTime();

        // Execute a Runge Kutta step
        auto step = runge_kutta.step();

        // Get the current result and timestep
        auto timestep = runge_kutta.getTimeStep();
        position = runge_kutta.getValue();

        // Get electric field at current position and fall back to empty field if it does not exist
        auto efield = detector_->getElectricField(static_cast<ROOT::Math::XYZPoint>(position));

        // Apply diffusion step
        auto diffusion = carrier_diffusion(std::sqrt(efield.Mag2()), timestep);
        position += diffusion;
        runge_kutta.setValue(position);

        // Adapt step size to match target precision
        double uncertainty = step.error.norm();

        // Lower timestep when reaching the sensor edge
        if(std::fabs(model_->getSensorSize().z() / 2.0 - position.z()) < 2 * step.value.z()) {
            timestep *= 0.75;
        } else {
            if(uncertainty > target_spatial_precision_) {
                timestep *= 0.75;
            } else if(2 * uncertainty < target_spatial_precision_) {
                timestep *= 1.5;
            }
        }
        // Limit the timestep to certain minimum and maximum step sizes
        if(timestep > timestep_max_) {
            timestep = timestep_max_;
        } else if(timestep < timestep_min_) {
            timestep = timestep_min_;
        }
        runge_kutta.setTimeStep(timestep);
    }

    // Find proper final position in the sensor
    auto time = runge_kutta.getTime();
    if(!detector_->isWithinSensor(static_cast<ROOT::Math::XYZPoint>(position))) {
        auto check_position = position;
        check_position.z() = last_position.z();
        if(true) { // position.z() > 0 && detector_->isWithinSensor(static_cast<ROOT::Math::XYZPoint>(check_position))) {
            // Carrier left sensor on the side of the pixel grid, interpolate end point on surface
            auto z_cur_border = std::fabs(position.z() - model_->getSensorSize().z() / 2.0);
            auto z_last_border = std::fabs(model_->getSensorSize().z() / 2.0 - last_position.z());
            auto z_total = z_cur_border + z_last_border;
            position = (z_last_border / z_total) * position + (z_cur_border / z_total) * last_position;

            time = (z_last_border / z_total) * time + (z_cur_border / z_total) * last_time;

            // Check if position is outside the detector
            if(std::fabs(position.z()) > model_->getSensorSize().z() / 2.0) {
                position.z() = model_->getSensorSize().z() / 2.0;
                if(position[2] < 0) {
                    position[2] *= -1;
                }
            }
        } else {
            // Carrier left sensor on any order border, use last position inside instead
            position = last_position;
            time = last_time;
        }
    }

    double total_time = time; // runge_kutta.getTime();

    // Apply repulsion to the charges
    double epsilon_0r = 11.68 * 8.854187817e-12 * Units::get(1., "m");
    double repulsion_coeff = carrier_mobility_init*charge_density_ / (ROOT::Math::Pi()*epsilon_0r) * total_time;

    double repulsion_std_dev;
    if(cylinder) {
        repulsion_std_dev = std::sqrt(repulsion_coeff);
    } else {
        repulsion_std_dev = std::pow(repulsion_coeff, 2./3);
    }

    // Compute the repulsion in three dimensions
    Eigen::Vector3d repulsion;
    std::uniform_real_distribution<double> radius_distribution(0, repulsion_std_dev);
    std::uniform_real_distribution<double> phi_distribution(0., 2*ROOT::Math::Pi());

    // cylinder = true;
    if(cylinder) {
        std::uniform_real_distribution<double> radius_squared_distribution(0, repulsion_std_dev*repulsion_std_dev);
        double rand_rep = radius_squared_distribution(random_generator_);
        double phi = phi_distribution(random_generator_);

        // Get the minimum component
        auto abs_v_vec = v_vec;
        for(auto& elm : v_vec) {
           elm = std::abs(elm); 
        }
        auto minIt = std::min_element(abs_v_vec.begin(), abs_v_vec.end());
        int idxMin = static_cast<int>(minIt - abs_v_vec.begin());

        std::rotate(v_vec.begin(), v_vec.begin()+idxMin, v_vec.end());

        // Initialize placeholders
        double A = v_vec.at(0);
        double B = v_vec.at(1);
        double C = v_vec.at(2);

        // Orthonormal vectors
        Eigen::Vector3d u = Eigen::Vector3d(0., -C, B);
        Eigen::Vector3d v = Eigen::Vector3d(B*B + C*C, -A*B, -A*C);
        u /= u.norm();
        v /= v.norm();

        // LOG(INFO) << u << " " << std::cos(phi);
        // LOG(INFO) << v << " " << std::sin(phi);
        // https://stackoverflow.com/questions/5837572/generate-a-random-point-within-a-circle-uniformly

        repulsion = std::sqrt(rand_rep) * std::cos(phi) * u + std::sqrt(rand_rep) * std::sin(phi) * v;
        // repulsion = rand_rep * (u * std::cos(phi) + v * std::sin(phi));

        // Reset z-value
        repulsion[2] = 0;

    } else {
        std::uniform_real_distribution<double> theta_distribution(0., 1);

        double phi = phi_distribution(random_generator_);
        double theta = theta_distribution(random_generator_);
        double cosTheta = 1 - 2 * theta;
        double sinTheta = std::sqrt(1 - cosTheta*cosTheta);

        double rad = radius_distribution(random_generator_);
        repulsion = Eigen::Vector3d(rad * sinTheta*std::cos(phi), rad * sinTheta*std::sin(phi), rad * cosTheta);
        // LOG(INFO) << repulsion;
    }

    // LOG(INFO) << cylinder << " " << repulsion_std_dev;
    // LOG(INFO) << position;
    position += repulsion;
    // LOG(INFO) << position;
    // LOG(INFO) << "";

    // Return the final position of the propagated charge
    return std::make_pair(static_cast<ROOT::Math::XYZPoint>(position), time);
}

void GenericPropagationModule::finalize() {
    long double average_time = total_time_ / std::max(1u, total_propagated_charges_);
    LOG(INFO) << "Propagated total of " << total_propagated_charges_ << " charges in " << total_steps_
              << " steps in average time of " << Units::display(average_time, "ns");

    if(output_plots_) {
        drift_time_histo->Write();
    }
}
