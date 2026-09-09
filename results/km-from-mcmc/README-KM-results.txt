Subgroups use the fitted posterior-mean gamma values and the package score definition: 1 + biomarkers %*% gamma.
Scores greater than or equal to zero are classified as above threshold.
Kaplan-Meier curves and log-rank tests describe survival differences but do not directly test beta 3.
The beta-3 Wald and likelihood-ratio p-values come from Cox models with treatment, subgroup, and treatment-by-subgroup interaction.
Because gamma and subgroup membership were estimated from these same observations, all p-values are exploratory/post-selection and may be anti-conservative.
A confirmatory analysis requires a prespecified score evaluated in independent validation data or a method that propagates subgroup-estimation uncertainty.
