.. |bullet| unicode:: U+02022
.. |emdash| unicode:: U+02014

========================================
Classification for Renal Cell Carcinomas
========================================

Final project presentation
~~~~~~~~~~~~~~~~~~~~~~~~~~

Goal
====

Classification of

  * Kidney renal papillary cell carcinoma (KIRP)
  * Kidney renal clear cell carcinoma (KIRC)

based on TCGA RNAseq abundance estimates


Data: Expected read counts with RSEM
====================================

.. image:: ../../fig/rsem.png
   :align: center

Data: training and validation sets
==================================

* 72 samples

  * 35 KIRC / 37 KIRP
  * 47 training / 25 test
  * 24 KIRC / 23 KIRP (training)

Method
======

EDA

* Log-transformed
* Filtering
* Normalization

Feature Selection

* Variance
* T-score

Classification methods

* LDA with features selected
* SVM with features selected
* Random forest


EDA: Filtering
==============

.. image:: ../../fig/filter_histograms.png
   :align: center

EDA: Normalization
==================

.. image:: ../../fig/postfilter_boxplot.png
   :align: center

EDA: Normalization
==================

.. image:: ../../fig/postnorm_test_boxplot.png
   :align: center

EDA: Normalization
==================

.. image:: ../../fig/postnorm_histogram.png
   :align: center

EDA: MD plots
=============

.. image:: ../../fig/mdplots_postnorm.png
   :align: center

EDA: PCA
========

.. image:: ../../fig/pcplot.png
   :align: center

Method
======

Feature Selection

* Variance
* T-score

Classification methods

* LDA with features selected
* SVM with features selected
* Random forest

Results: Performance
====================

.. image:: ../../fig/testset_acc.png
   :align: center

Results: LDA confusion matrix
=============================

::

     Variance                 Tscore
    
        true                    true
  pred   KIRC KIRP        pred   KIRC KIRP
    KIRC    7    4          KIRC    8    3
    KIRP    4   10          KIRP    3   11

Results: SVM confusion matrix
=============================

::

     Variance                 Tscore
    
        true                    true
  pred   KIRC KIRP        pred   KIRC KIRP
    KIRC    7    0          KIRC   10    0
    KIRP    4   14          KIRP    1   14

Results: Random Forest confusion matrix
=======================================

::

                   true
             pred   KIRC KIRP
               KIRC    8    0
               KIRP    3   14 

Summary
=======

* Good performance
* Test statistic worked better than variance ranking
* Future directions

  * Technical and/or biological nuissance effects
  * Error bars
  * Biological significance
  * Sensitivity study
  * Non-cancer data
  
