clc;clear;
codeRates = ["1-2";"1-3";"22-26";"22-27"];
foldername ="LDPC/Minsum/";
for cc = 1:numel(codeRates)
    cfoldername = foldername+"CR"+codeRates(cc);
    d = dir(cfoldername);
    fnames = {d(3:end).name};
  
    for kk =1:numel(fnames)
       t{kk}.datatable=readtable( cfoldername+"/"+fnames{kk});
       t{kk}.codeRates=codeRates(cc);
       tmp = strsplit(fnames{kk},{'-','.'});
       Q = str2double(tmp{end-2}(2:end));
       M =  str2double(tmp{end-1}(2:end));
       t{kk}.Q = Q;
       t{kk}.M = M;
    end
end
