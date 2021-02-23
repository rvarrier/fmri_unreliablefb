cvMANOVA = load('cvMANOVA_V1results.mat');
V1_PD = cvMANOVA.V1_PD;
% convert to table
d = array2table(V1_PD.data,'VariableNames',[V1_PD.labels(1:4)' {'Fbtype'},{'deviation'}]);
d.subject = (1:size(d,1))';
d.subject(d.subject>21) = d.subject(d.subject>21)+1; % subject 22 was removed;
d.thr = thr(d.subject);
d = stack(d,{'Pre_int_stimpair_1'    'Pre_int_stimpair_2'    'Post_int_stimpair_1'    'Post_int_stimpair_2'},'NewDataVariableName','cv','IndexVariableName','condition');

d.phase = cellfun(@(x)any(regexpi(x,'post.*')),cellstr(d.condition));
label = {'1_pre','3_post'};
d.phase = label(d.phase+1)';
d.stim = cellfun(@(x)any(regexpi(x,'_2')),cellstr(d.condition))+1;

fitlme(d,'cv~1+phase*Fbtype + (1+phase|subject)')

fitlme(d,'cv~1+thr+phase*Fbtype + (1+phase|subject)')
% Save it as easily accessible csv
writetable(d,'d_cv.csv')

%% plotting using gramm-toolbox
g = gramm('y',d.cv,'x',d.phase,'color',d.Fbtype); 
g.geom_point('alpha',0.1);

g.geom_line('alpha',0.1)
g.stat_summary('type','bootci')
g.draw()


%
figure
g = gramm('y',d.cv,'x',d.phase,'color',d.Fbtype); 
g.geom_point('alpha',0.1);

g.geom_line('alpha',0.1)
g.stat_summary('type','bootci')
g.facet_grid(d.stim,[])
g.draw()
%% Mixed Model
fitlme(d,'cv~-1+phase*Fbtype + (-1+phase|subject)')

fitlme(d,'cv~stim*phase*Fbtype + (1+stim*phase|subject)')

%%
tmp = load('data_across_blocks_and_subjects.mat');
% add FbType & subjectnumber
for sub = 1:size(tmp.Aalldata,1)
    for run = 1:size(tmp.Aalldata,2)
        old = tmp.Aalldata{sub,run};
        tmp.Aalldata{sub,run} =[ old ...
            repmat(sub,size(old,1),1) ...
            repmat(tmp.Afbtype(sub),size(old,1),1) ...
            repmat(tmp.thr(sub),size(old,1),1)];
    end
end
% labellists
% there are some undescribed unknowns
preList = {'run','trial','stim','response','response_unmapped','trialOnset','stimOnset','respOnset','subject','Fbtype','threshold'};
feedList = {'run','trial','stim','response','response_unmapped','possibly_staircase','unknown1','unknown2','unknown_timing1','unknown_timing2','unknown_timing3','unknown_timing4','subject','Fbtype','threshold'};

% generate 3 data frames for the respective phases
%pre
x = tmp.Aalldata(:,1:8);
pre = array2table(cat(1,x{:}),'VariableNames',preList);
pre.phase = repmat({'1_pre'},1,size(pre.run,1))';
%feed
x = tmp.Aalldata(:,9:16);

feedlist= array2table(cat(1,x{:}),'VariableNames',feedList);
feedlist.run = feedlist.run + 8;
feedlist.phase = repmat({'2_feedback'},1,size(feedlist.run,1))';
%post
x = tmp.Aalldata(:,17:end);
post= array2table(cat(1,x{:}),'VariableNames',preList);
post.run = post.run + 16;
post.phase = repmat({'3_post'},1,size(post.run,1))';
% concatenate them, remove unknown & unmatchable fields
d_beh = [pre(:,[1:5 end-3:end]); feedlist(:,[1:5 end-3:end]); post(:,[1:5 end-3:end])]; 

% This I think should be correct, copied from behavioral_performance.m
d_beh.correct = mod(d_beh.stim,2) == (d_beh.response==1);

% Again I think this should be correct
label = {'reliable','unreliable'};
d_beh.Fbtype_str  = label(d_beh.Fbtype)';

% Save it as easily accessible csv
writetable(d_beh,'d_beh.csv')
%% Some plotting
g = gramm('y',d_beh.correct,'x',d_beh.phase,'color',d_beh.Fbtype_str)
g.stat_summary('geom','errorbar','type','bootci')
g.draw()
%% run  mixed effects GLM
d_beh.stim_cat = categorical(d_beh.stim);
fitglme(d_beh,'correct~1+stim_cat+Fbtype+(1+stim_cat|subject)','Distribution','Binomial')


