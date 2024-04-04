# Multimodal Biomedical Data Fusion Using Sparse Canonical Correlation Analysis and Cooperative Learning: A Cohort Study on COVID-19

Ahmet Gorkem Er<sup>1,2,3</sup>, Daisy Yi Ding<sup>4</sup>, Berrin Er<sup>5</sup>, Mertcan Uzun<sup>3</sup>, Mehmet Cakmak<sup>6</sup>, Christoph Sadee<sup>1</sup>, Gamze Durhan<sup>7</sup>, Mustafa Nasuh Ozmen<sup>7</sup>, Mine Durusu Tanriover<sup>6</sup>, Arzu Topeli<sup>5</sup>, Yesim Aydin Son<sup>2</sup>, Robert Tibshirani<sup>4,8</sup>, Serhat Unal<sup>3</sup>, Olivier Gevaert<sup>1,4</sup>

1)	Stanford Center for Biomedical Informatics Research (BMIR), Department of Medicine, Stanford University, Stanford, CA, 94305, USA
2)	Department of Health Informatics, Graduate School of Informatics, Middle East Technical University, Ankara, 06800, Türkiye
3)	Department of Infectious Diseases and Clinical Microbiology, Hacettepe University Faculty of Medicine, Ankara, 06230, Türkiye
4)	Department of Biomedical Data Science, Stanford University, Stanford, CA, 94305, USA
5)	Department of Internal Medicine, Division of Intensive Care Medicine, Hacettepe University Faculty of Medicine, Ankara, 06230, Türkiye
6)	Department of Internal Medicine, Hacettepe University Faculty of Medicine, Ankara, 06230, Türkiye
7)	Department of Radiology, Hacettepe University Faculty of Medicine, Ankara, 06230, Türkiye
8)	Department of Statistics, Stanford University, Stanford, CA, 94305, USA

Corresponding authors:

Ahmet Gorkem Er
E-mail: ahmetgorkemer@gmail.com

Olivier Gevaert
E-mail: ogevaert@stanford.edu

## Abstract

Through technological innovations, patient cohorts can be examined from multiple views with high-dimensional, multiscale biomedical data to classify clinical phenotypes and predict outcomes. Here, we aim to present our approach for analyzing multimodal data using unsupervised and supervised sparse linear methods in a COVID-19 patient cohort. This prospective cohort study of 149 adult patients was conducted in a tertiary care academic center. First, we used sparse canonical correlation analysis (CCA) to identify and quantify relationships across different data modalities, including viral genome sequencing, imaging, clinical data, and laboratory results. Then, we used cooperative learning to predict the clinical outcome of COVID-19 patients: ICU admission. We show that serum biomarkers representing severe disease and acute phase response correlate with original and wavelet radiomics features in the LLL frequency channel (𝑐𝑜𝑟𝑟(𝑋u𝟏, Zv𝟏) = 0.596, p-value < 0.001). Among radiomics features, histogram-based first-order features reporting the skewness, kurtosis, and uniformity have the lowest negative, whereas entropy-related features have the highest positive coefficients. Moreover, unsupervised analysis of clinical data and laboratory results gives insights into distinct clinical phenotypes. Leveraging the availability of global viral genome databases, we demonstrate that the Word2Vec natural language processing model can be used for viral genome encoding. It not only separates major SARS-CoV-2 variants but also allows the preservation of phylogenetic relationships among them. Our quadruple model using Word2Vec encoding achieves better prediction results in the supervised task. The model yields area under the curve (AUC) and accuracy values of 0.87 and 0.77, respectively. Our study illustrates that sparse CCA analysis and cooperative learning are powerful techniques for handling high-dimensional, multimodal data to investigate multivariate associations in unsupervised and supervised tasks.

