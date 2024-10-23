
if exist(exp_path)==7
    inp = input('this expirement exists. Delete it? (y/n): ','s');
    if strcmp(inp,'y')
        cmd_rmdir(exp_path);
    elseif strcmp(inp,'n')
        fprintf('Not deleting')
        return
    else
        fprintf('Wrong input\n')
        error('this expirement exits.')
    end
end
mkdir(exp_path);
mkdir([exp_path '/results']);
mkdir([exp_path '/figs']);

fprintf(['Check ' exp_path '/expirement_info.txt for more info\n'])

fid = fopen([exp_path '/expirement_info.txt'],'wt');
info='This expirement will produce RS and vanilla SISO system figures\n';
info = [info '\n---System parameters--- \n'];
info = [info 'Results export Path: ' exp_path '\n'];
info = [info 'cnq path export Path: ' cnq_path '\n'];
info = [info 'Run Parallel: ' char(string(run_parallel)) '\n'];

info = [info '\n---Experiment parameters--- \n'];
Qrsmat = Qmat;
SS=size(PbitB,2);
Pbit_lin = 10.^(PbitB/10);
info = [info 'Modulation orders: ' num2str(Qmat') '\n'];
info = [info 'Noise Modes: ' num2str(allM') '\n'];
info = [info 'Pbit =[' num2str(PbitB(1)) ',' num2str(PbitB(end)) '] with step ' num2str(step) '\n'];


%RS(255,239)
info =[info '\n---RS params ----\n'];
info = [info 'Number of RS batches: ' num2str(N) '\n'];
info =[info 'n=' num2str(nn_pair') '\n'];
info =[info 'k=' num2str(kk_pair') '\n'];
mm_pair  = log2(nn_pair+1);
info =[info 'm=' num2str(mm_pair') '\n'];
info =[info 'd=' num2str(nn_pair'-kk_pair' +1) '\n'];
info =[info 't=' num2str(floor((nn_pair'-kk_pair')./2)) '\n'];
info=[info 'Use Interleving =' char(string(useInterleving)) '\n'];
fprintf(fid, info);
fclose(fid);

berCell = cell(length(Qmat),length(allM),length(kk_pair));
count = 0;

% if run_parallel && exist('max_cores', 'var') && ~isempty(max_cores)
%         parpool('local')
% end