function genCCode(f, fName, args)
% generate C code from matlab functions

% input:
%  f: symbolic Function
% fName: name of the function
% args: a set of objects whose length are used to compute the length of the
% arguments

N = numel(f); % size of output argument

% 	argsList = ', '.join([('const double *' if isinstance(q,tuple) else 'const double ') + p for p,q in args])
argsList = [];

% input arguments
for i=1:length(args)
    if length(args{i}) > 1
        % a vector
        s = sprintf('const double in%d[%d], ', i, length(args{i}));
    else
        % a scalar
        s = sprintf('const double in%d, ', i);

    end
    argsList = [argsList, s];
end

% generate the header
header = [sprintf('#ifndef %s_H\n', upper(fName)), ...
    sprintf('#define %s_H\n', upper(fName)), ...
    sprintf('%s\n','#include <math.h>'), ...
    sprintf('%s\n','#include <stdlib.h>'), ...
    sprintf('%s\n','#include <stdio.h>'), ...
    sprintf('%s\n','#include <string.h>'), ...
    sprintf('void %s(%s double varOut[%d]);\n', fName, argsList, N), ...
    '#endif'];

% disp(header)

writeText(header, sprintf('../C/%s.h', fName))

% generate code

f=f.'; % use transpose of the symbolic object so that it's in the usual matlab column-major ordering

code = [sprintf('#include \"%s.h\"\n', fName), ...
    sprintf('void %s(%s double varOut[%d])\n{\n', fName, argsList, N), ...
    sprintf('double (*A0)[%d][%d] = varOut;\n', size(f,1), size(f,2)), ...
    sprintf('memset(A0, 0.0, sizeof(double)*%d);\n', N)];

% c version of the function
% c = ccode(f); 

ccode(f, 'file', sprintf('../C/%s_tmp.c', fName));

lines = fileread(sprintf('../C/%s_tmp.c', fName));
lines = regexprep(lines, '(t\d+ =)', 'double $1'); % declare temps
lines = strrep(lines, 'A0', '(*A0)'); % declare temps
delete(sprintf('../C/%s_tmp.c', fName));
% varname = @(x) inputname(1);

% get the variable names
for i=1:numel(args)
    if length(args{i})==1
        continue
    end
    
    for j=1:length(args{i})
        code = [code, sprintf('double %s = in%d[%d];\n', char(args{i}(j)), i, j-1)];
    end
end

code = [code, lines];


% copy the result to the output variable
code = [code, '}']; %sprintf('\nmemcpy(varOut, A0, sizeof(double)*%d);\n', N), '}'];

% disp(code)

writeText(code, sprintf('../C/%s.c', fName))

% def genCcode(f, fName, args ):
% 	'''generate a very simple C program from the matrix function passed in and the
% 	arguments necessary'''
% 
% 	argsList = ', '.join([('const double *' if isinstance(q,tuple) else 'const double ') + p for p,q in args])
% 	N = len(f)	
% 
% 	"header file"
% 	header = """
% #ifndef %s_H
% #define %s_H
% #include "math.h"	
% void %s(%s, double varOut[%d]);
% #endif""" % (fName.upper(), fName.upper(), fName, argsList, N)
% 	
% 	#print header
% 	"source code"
% 	code = "#include \"%s.h\"\n" % fName
% 	code += "void %s(%s, double varOut[%d])\n" % (fName, argsList, N)
% 	#x = sum(f.reshape(len(f),1).tolist(),[])
% 	print(f)
% 	x = sum(f.transpose().reshape(len(f),1).tolist(),[])
% 
% 	code += "{\n"
% 	for name, value in args:
% 		k = 0
% 		"unpack local variables"
% 		if isinstance(value,tuple):
% 			for v in value:			
% 				code += "double %s = %s[%d];\n" % (v, name, k)
% 				k += 1
% 			
% 	for i in range(N):
% 		print(x[i])
% 		#xOut = parseString(str(x[i]))				
% 		code += "varOut[%d] = %s;\n" % (i, printing.ccode(x[i]))
% 		
% 	code += "}"
% 	#print code
% 	
% 	return (header, code)