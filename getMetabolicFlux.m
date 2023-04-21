%% To get metabolic fluxes from expresison data/raw TPM/FPKM/Counts

% upload human model and find core rxns (lactate rxns added later on)
load('Recon3E_O2Fixed.mat') % human model
model = modelO2Fixed
[grRatio, grRateKO, grRateWT, hasEffect, delRxn, fluxSolution] = singleRxnDeletion(model,90);
core = find(grRatio<0.5);

% Import exp data file
expData= readtable(**File Name**);
genes= cellstr(string(expData.Var1)) %Var1 when column has no header
expData.Var1=[];
names= expData.Properties.VariableNames;

% get and save models using fastcoreWeighted_mod
for i=1:width(expData)
  levels = table2array(expData(:,i));
  reaction_levels = gene_to_reaction_levels(model,genes, levels , @min, @(x,y)(x+y));
  reaction_levels(isnan(reaction_levels))=0;
  reactionWeights = max(reaction_levels)-reaction_levels;
  tissueModel = fastCoreWeighted_mod(core, model, reactionWeights);	
  savename = names{i}
  save(savename, 'tissueModel')
end

%  loop loading of models 
for i =1:width(expData)
  modelName = names{i}
  matfilename = strcat(modelName,'.mat');
  model = importdata(matfilename);

  % constrain biomass 10% and remove objective
  sol = optimizeCbModel(model)
  ind = find(ismember(model.rxns,'biomass_reaction'))
  model.lb(ind) = 0.1*sol.f
  model.c(model.c==1)=0;

  %Preprocessing
  modelEP = pre_processing(model);

  % Set sampling parameters
  exp_i = 0;
  av_exp = 0;
  va_exp = 0;
  Beta=1e8;
  damping=0.9;
  precision=1e-5;
  minvar=1e-50;
  maxvar=1e50;
  maxit=1e3;

  % run metabolicEP to get flux
  [mu_free_un, s_free_un, a_free_un, d_free_un, av_free_un, va_free_un, cov_free_un, t_EP_free_un] = ...
  MetabolicEP(full(modelEP.S), modelEP.b, modelEP.lb, modelEP.ub, Beta, damping, maxit, minvar, maxvar, precision, av_exp, va_exp, exp_i);
  modelind = find(ismember(modelO2Fixed.rxns,model.rxns));
  commetnew = cell(6465,1);
  commetnew(modelind)=num2cell(av_free_un);
  newMatrix(:,i) = commetnew;
end

% save flux into csv
newMatrix(cellfun('isempty',newMatrix))={0};
newTable = cell2table(newMatrix);
newTable Properties.VariableNames = names;
newTable.Properties.RowNames = modelO2Fixed.rxns;
writetable(newTable,'flux.csv','WriteRowNames',true);
