# RFES
Recursive Feature Elimination with Strikeout


The RFES algorithm proceeds as follows:

s and Accuracymax are initialised to 0.

The feature importance for each gene across all Groups in datatrain is obtained using the varImp function from the caret R package, and the maximum variable importance, as Mean Decrease in Accuracy (MDA) for each gene (j) across the Groups (m) is recorded. 
〖MDA〗_(j,max)=max⁡(〖MDA〗_(j,1),〖MDA〗_(j,2),…,〖MDA〗_(j,m) )

The first quantile (Q1) is found from the maximum VI across all features, where p is the total number of features. 
Q_(1  )= 〖Quantile〗_0.25 (〖MDA〗_(1,max),〖MDA〗_(2,max),…,〖MDA〗_(p,max))

The features (F) are retained if they have a maximum MDA greater than the first quantile.
F=(j ┤|  〖MDA〗_(j,max)>Q_1)

Using this feature set F, the model is run again and the accuracy of the model on datatrain and datatest are recorded as follows:
Accuracy=  (True Positives+True Negatives)/(True Positives+True Negatives+False Positives+False Negatives)

If Accuracytest (model accuracy on datatest) is greater than Accuracy max then Accuracy max =  Accuracy test, else s is incremented by 1. This is repeated until s > strikeout.

Repeat the procedure rpt times and record the most accurate model for each repeat.

Published in _Unsupervised machine learning identifies early stage breast cancer subclasses from aberrant glycogene expression and informs personalised drug treatment_ (https://doi.org/10.21203/rs.3.rs-10155568/v1)

