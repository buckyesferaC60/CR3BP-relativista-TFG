function Orbit = Read_Orbit_final(Orbit_struct)
%READ_ORBIT Orbit_struct:
%   .mu
%   .kappa
%   .family
%   .ci

    if Orbit_struct.family == 'Ly_L3'
        decimal = 5;
    else
        decimal = 6;
    end

    Orbit_dir = strrep(strcat('Orbit_files_final','\',num2str(Orbit_struct.kappa,10),'\',num2str(Orbit_struct.mu,10),'\',Orbit_struct.family),'.','p');

    Previous_dir = cd(Orbit_dir);

    %Orbit_filename = strcat(strrep(strcat('Orb',num2str(Orbit_struct.ci,20)),'.','p'),'.txt');
    
    if ~isa(Orbit_struct.ci,'double')
        Orbit_filename = strcat('Orb',Orbit_struct.ci,'.txt');
    else
        if abs(Orbit_struct.ci) < abs(round(Orbit_struct.ci,decimal))
            Orbit_struct.ci = Orbit_struct.ci - 10^-(decimal);
        end
        

        if isfield(Orbit_struct,'ci2')
            if abs(Orbit_struct.ci2) < abs(round(Orbit_struct.ci2,decimal))
                Orbit_struct.ci2 = Orbit_struct.ci2 - 10^-(decimal);
            end
            Orbit_filename = strcat(strrep(strcat('*','&',num2str(Orbit_struct.ci2,decimal+1)),'.','p'),'*','.txt');
        else
            Orbit_filename = strcat(strrep(strcat('Orb',num2str(Orbit_struct.ci,decimal+1)),'.','p'),'*','.txt');
        end   
    end

    Orbit_filename_struct = dir(Orbit_filename);

    Orbit = readmatrix(Orbit_filename_struct(1).name);

    cd(Previous_dir);


end
