function [subject, session, task, run] = acquisition_names(filename)

subject = [];
session = [];
task = [];
run = [];

strfile = strsplit(filename, {'_', '.'});
disp(['>>>>>>> Reading file: ', filename]);
for k = 1:length(strfile)
    
    detail = strfile{k};
    strdetail = strsplit(detail, {'-'});
    detail_type = strdetail{1};
    
    switch detail_type 
        
        case 'task'
            task = detail;
            
        case 'sub'
            subject = detail;
            
        case 'ses'
            session = detail;
        
        case 'run'
            run = detail;
        
        case 'eeg'
            continue
            
        otherwise 
            disp('Unrecognized field, please check the name of your file')
    end
end