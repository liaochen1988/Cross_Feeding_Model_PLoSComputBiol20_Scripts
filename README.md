# A coarse-grained ecology model for microbial communities

Our paper has been published in PLoS Computational Biology (https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008135).

The project aims to advance mathematical modeling of microbial communities towards applications to realistic data. Most of the existing models either lack mechanistic details or too complex to understand. We propose that an intermediate-scale model may be a better solution to infer the underlying microbial processes that drive temporal changes in microbial community state (e.g., population density, nutrient concentration, ...). We developed such a coarse-grained ecology model by combining population dynamics models and a simplifed metabolite network for each cell type (or ecotype) in the community. This repository includes our MATLAB codes for three cross-feeding examples that demonstrate the utility of our new framework in filling knowledge gaps in microbial ecology. The three cross-feeding ecosystems include

(1) Unilateral acetate-mediated cross-feeding between a glucose specialist and an acetate specialist. Ref: Rosenzweig, R. Frank, et al. "Microbial evolution in a simple unstructured environment: genetic differentiation in Escherichia coli." Genetics 137.4 (1994): 903-917.

(2) Bilateral amino-acid-mediated cross-feeding between a lysine auxotroph and a leucine auxotroph. Ref: Zhang, Xiaolin, and Jennifer L. Reed. "Adaptive evolution of synthetic cooperating communities improves growth performance." PloS one 9.10 (2014): e108297.

(3) Multilateral amino-acid-mediated cross-feeding between 14 amino acid auxotrophs. Ref: Mee, Michael T., et al. "Syntrophic exchange in synthetic microbial communities." Proceedings of the National Academy of Sciences 111.20 (2014): E2149-E2156. 
