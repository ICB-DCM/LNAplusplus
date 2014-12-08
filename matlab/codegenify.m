% helper function to add %codegen tags to matlab-generated functions
% in order for the Coder toolbox to process them

function codegenify(srcDir)

% find .m files
files = dir([srcDir '/*.m']);

N = length(files);
for i=1:N
    f = fopen([srcDir '/' files(i).name], 'r');
    f2 = fopen('tmp','w');
    
    while ~feof(f)
        l = fgetl(f);
        if strmatch('function ', l)
            if isempty(strfind(l, '%#codegen'))
                l = [l ' %#codegen'];
            end
        end
        fprintf(f2, '%s\n', l);
%         fprintf('%s\n', l);
    end
    fclose(f);
    fclose(f2);
    movefile('tmp', [srcDir '/' files(i).name]);
end

