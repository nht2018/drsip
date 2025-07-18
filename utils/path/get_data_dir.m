
function data_dir = get_data_dir(datasetname)
    if nargin < 1
        datasetname = '';
    else
        datasetname = char(datasetname);
    end
    if strcmp(datasetname, 'qplib')
         if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
            data_dir = '/Users/niehantao/Desktop/software/data/sdp_data/qplib/html/lp';
        elseif strcmp(computer, 'GLNXA64') % my linux server
           error('not implemented yet');
        end
    elseif strcmp(datasetname, 'NETLIB_MAT_MDO_PRSLVD') || strcmp(datasetname, 'MITTELMANN_MAT_MDO_PRSLVD')
        if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
            data_dir = '/Users/niehantao/Desktop/SSNLP/data';
        elseif strcmp(computer, 'GLNXA64') % my linux server
           error('not implemented yet');
        end
    else
        if strcmp(computer, 'MACI64') || strcmp(computer, 'MACA64') % my mac laptop
            data_dir = ['/Users/niehantao/Desktop/software/data/sdp_data/', datasetname, '/'];
        elseif strcmp(computer, 'GLNXA64') % my linux server
            data_dir = ['/public/shared/sdp_data/', datasetname, '/'];
        end
    end
end
