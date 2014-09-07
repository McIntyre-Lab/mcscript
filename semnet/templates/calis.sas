libname SEM "$libname";

proc calis data=SEM.$data method=fiml COVARIANCE RESIDUAL MODIFICATION maxiter=10000 outfit=fitstat;
    lismod
        yvar = $yvar, 
        xvar = $xvar
        ;
        matrix _BETA_ $beta
        ;
        matrix _GAMMA_ $gamma
        ;
        matrix _PHI_ $phi
        ;
    run;

data SEM.gene_${gene}_model_${model};
    length gene $$12.;
    length model $$14.;
    retain gene model;
    set fitstat;
    where IndexCode = 312;
    gene = "$gene";
    model = "Model $model";
    rename FitValue = BIC;
    keep gene model FitValue;
    run;
