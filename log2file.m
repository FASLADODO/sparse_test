%this script is used to log console output into a log file
function log2file(fid, content)
%content, the log content to be stored
if fid
    fprintf(fid, '%s\n', content);
end

end