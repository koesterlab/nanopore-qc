Plot of k-mer hits obtained with Kraken_ for the first 100 bases of the first 100 reads.
Each pixel represents a particular k-mer in one read, which is colored by the taxon it maps to.
Taxon colors are the same as in the classification tree.
They are obtained by

1. Calculating a circular layout for the unfiltered classification tree.
2. Mapping angle and magnitude of each node to HSV_ color space.

This procedure ensures that

1. Taxons yield a more saturated color the further away (i.e. more specific) they are from the root node.
2. Taxons on the same branch get a similar hue.

Unclassified k-mers or k-mers mapping to taxa that do not have enough abundance to occur in the classification tree are shown as black.

.. _Kraken: http://ccb.jhu.edu/software/kraken
.. _HSV: https://en.wikipedia.org/wiki/HSL_and_HSV

