function filename = Generate_Orbit_file(str,Orbit)
%GENERATE_ORBIT_FILE Summary of this function goes here
%   Detailed explanation goes here

    mu_string = num2str(str.mu);
    kappa_string = num2str(str.kappa);
    x0_string = num2str(str.x0);
    z0_string = num2str(str.z0);

    filename = strrep(strcat(...
        'mu',mu_string,'_',...
        'K',kappa_string,'_',...
        'Type',str.family_type,str.family_N,'_',...
        'X',x0_string,'_',...
        'Z',z0_string),'.','p');

    filename = strcat(filename,'.txt');

    fileID = fopen(filename,'w'); % read r write w 

    writematrix(Orbit,filename);

    fclose(fileID);
    
end

