/*
 * Temporary producer of histogram hits
 */

#ifndef ALLPIX_MODULE_DETECTOR_HISTOGRAMMER_H
#define ALLPIX_MODULE_DETECTOR_HISTOGRAMMER_H

#include <memory>
#include <string>

#include <TFile.h>
#include <TH1I.h>
#include <TH2I.h>

#include "core/config/Configuration.hpp"
#include "core/geometry/GeometryManager.hpp"
#include "core/messenger/Messenger.hpp"
#include "core/module/Module.hpp"

#include "objects/PixelCharge.hpp"

namespace allpix {
    class DepositionMessage;

    // define the module to inherit from the module base class
    class DetectorHistogrammerModule : public Module {
    public:
        // constructor and destructor
        DetectorHistogrammerModule(Configuration, Messenger*, std::shared_ptr<Detector>);
        ~DetectorHistogrammerModule() override;

        // create the histograms
        void init() override;

        // fill the histograms
        void run(unsigned int) override;

        // write the histograms
        void finalize() override;

    private:
        // configuration for this module
        Configuration config_;

        // attached detector and model
        std::shared_ptr<Detector> detector_;

        // deposits for a specific detector
        std::shared_ptr<PixelChargeMessage> pixels_message_;

        // histograms
        TH2I* histogram; // FIXME: bad name
        TH1I* cluster_size;

        // output file
        TFile* output_file_;
    };
} // namespace allpix

#endif /* ALLPIX_MODULE_DETECTOR_HISTOGRAMMER_H */