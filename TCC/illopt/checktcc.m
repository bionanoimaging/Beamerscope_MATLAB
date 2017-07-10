function [ filepath_dbtcc tcci tccavailable] = checktcc(IOparams, database_params)
%checktcc searches for TCCs which have been precalculated before
%%% input
%   IOparams - set of parameters for current optimization task
%   database_params - database with parameters which have been
%   precalculated already
%%% output
%   filepath_dbtcc - filepath to the TCC-database
%   tcci - iterator of the database

for ( i=1:size(database_params,2) )
    
    % check if parameters have been used before
    if(isequal(IOparams, database_params(i)))
        tcci = i;
        filepath_dbtcc = strcat(database_params(i).path_db_tcc, num2str(tcci));
        display('Dataset has been found!')
        tccavailable = true; 
        return
    end
    
end

tcci = i + 1;
tccavailable = false;
filepath_dbtcc = strcat(IOparams.path_db_tcc, num2str(tcci));
display('Dataset has not been found!')

