function Out = Write_Orbit_final(Orbit_struct,Orbit)
%WRITE_ORBIT Orbit_struct:
%   .mu
%   .kappa
%   .family
%   .ci

    Orbit_dir = strrep(strcat('Orbit_files_final','\',num2str(Orbit_struct.kappa,10),'\',num2str(Orbit_struct.mu,10),'\',Orbit_struct.family),'.','p'); 

    if exist(Orbit_dir,"dir") ~= 7
        mkdir(Orbit_dir)
    end

    Previous_dir = cd(Orbit_dir);

    Orbit_filename = strrep(strcat('Orb',num2str(Orbit_struct.ci,20)),'.','p');

    if isfield(Orbit_struct,'ci2')
        Orbit_filename = strrep(strcat(Orbit_filename,'&',num2str(Orbit_struct.ci2,20)),'.','p');
    end

    fileID = fopen(Orbit_filename,'w');

    writematrix(Orbit,Orbit_filename);

    fclose(fileID);

    cd C:\Matlab

    cd ..

    exp_startup;

    cd(Previous_dir);

    Out = 0;

end

