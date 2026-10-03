"""Diagnosis proxy audit (design assumption A7).

Diagnosis is absent from the targets but can leak through evidence features
(drug class, disease identifiers behind HPO annotations). The audit trains a
simple classifier to predict (a) ATC drug class and (b) OMIM disease identifier
from the evidence feature table and reports held-out AUROC. A high AUROC is a
finding to report and a reason to drop the leaking feature, not a modelling
choice to hide.

Not implemented in version 0.1.
"""


def audit_diagnosis_proxies(evidence_feature_table, drug_class_labels, disease_identifier_labels):
    raise NotImplementedError("fit logistic regression per label set with grouped folds; return AUROC per label set")
