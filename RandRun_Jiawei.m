tic
clear all

for zztop=1:1

lengthfactor = 0.54 
lff = lengthfactor;
LM = lengthfactor*[15.3408,16.2173,20.1895,20.8322,24.5774];
LS=0.2*LM;
LS(2)=1.2*LS(2);
LS(3)=0.9*LS(3);
LS(4)=0.8*LS(4);%LS =1.6568    2.1018    1.9624    1.7999    2.6544

l2_ptcdi = lengthfactor*25.5;
lmin_ptcdi = lengthfactor*17;
lp = lengthfactor*10.1;

b_alk = 1.47; %npArrayExperimental value      %/Angstrom alkane
b_ptc = 0.99; %npArrayExperimental value      %/Angstrom ptcdi

Vxmin = 18;
Vstep = 30;
Trials = 128;
NPgridrows = 16;
grids =1; %how many grids to generate

%% Initial Variables
% h1 is number of hexagon rows there are
h1 = NPgridrows;
H = h1*3^(1/2); % H is the hieght of the grid
W = 3*H; % W is the exact width needed to fulfill 3:1 ratio
w = round(W); % w is the width - using one bond length as the counting measure
b = w + 1; % b is the number of atoms in one of the short rows
c = b + 1; % c is the number of atoms in one of the long rows
na = b + (b + c)*h1; % na is the number of atoms
xm1 = 2*(w + 1) + 1; % xm1 is the maximum x position at the longer rows
xm2 = 2*(w + 1); % xm2 is the maximum x position at the shorter rows
ym = h1*2 + 1; % ym is the maximum y position
ymax = 2*ym - 1;

SmeG = zeros(grids,5);
HmeG = zeros(grids,5);
SmeGStd = zeros(grids,5);
HmeGStd = zeros(grids,5);
SmeGG = zeros(grids,5);
HmeGG = zeros(grids,5);
SmeGGStd = zeros(grids,5);
HmeGGStd = zeros(grids,5);
Matrix_OC = zeros(xm1, ymax, grids, 5);
Matrix_OCex = zeros(xm1, ymax, grids, 5);
Matrix_OLex = zeros(xm1, ymax, grids, 5);
Matrix_OL = zeros(xm1, ymax, grids, 5);
Matrix_OG = zeros(xm1, ymax, 2, grids, 5);
Matrix_OGex = zeros(xm1, ymax, 2, grids, 5);
Ea_alk=0.017;
Ea_ptc=0.0266;


Vx = 0.41


%% Main Loop

for ngrd=1:grids
    for j = 1:5
         [na, Ligand_Ex, Ligand, SmeGG(ngrd,j),HmeGG(ngrd,j), SmeGGStd(ngrd,j), HmeGGStd(ngrd,j), Matrix_Ex, G_Ex,SmeG(ngrd,j),HmeG(ngrd,j), SmeGStd(ngrd,j), HmeGStd(ngrd,j), Matrix, G] = Randl_Ex(NPgridrows, 298, LM(j), l2_ptcdi, lmin_ptcdi, LS(j), lp, b_alk, b_ptc, Ea_alk, Ea_ptc, Vx, Trials);


         for x = 1: (xm1)
                for y = 1: (ymax)
                     Matrix_OCex(x, y, ngrd, j) = Matrix_Ex(x, y); %Current
                     Matrix_OLex(x, y, ngrd, j) = Ligand_Ex(x, y); %Ligand Length
                     for z = 1:2
                        Matrix_OGex(x, y, z, ngrd, j) = G_Ex(x, y, z); %Conductance
                     end


                end
         end

         for x = 1: (xm1)
                for y = 1: (ymax)
                     Matrix_OC(x, y, ngrd, j) = Matrix(x, y);
                     Matrix_OL(x, y, ngrd, j) = Ligand(x, y);
                     for z = 1:2
                        Matrix_OG(x, y, z, ngrd, j) = G(x, y, z);
                     end


                end
         end

    end
end


HNGG = HmeGG./SmeGG;
HNG = HmeG./SmeG;

CNG = HNGG./HNG;
CNGs = std(HNGG./HNG);

timelapse = toc

%% Photocurrent & Ligand & Conductivity Matrix
Matrix_Cex = zeros(xm1, ymax, 5);
Matrix_Lex = zeros(xm1, ymax,5);
Matrix_C = zeros(xm1, ymax, 5);
Matrix_L = zeros(xm1, ymax, 5);
Matrix_G = zeros(xm1, ymax, 2, 5);
Matrix_Gex = zeros(xm1, ymax, 2, 5);

for j= 1:5
    for x = 1: (xm1)
        for y = 1: (ymax)
            Matrix_Cex(x, y, j) = mean(Matrix_OCex(x, y, :, j));
            Matrix_C(x, y, j) = mean(Matrix_OC(x, y, :, j));
            Matrix_Lex(x, y, j) = mean(Matrix_OLex(x, y, :, j));
            Matrix_L(x, y, j) = mean(Matrix_OL(x, y, :, j));

            for z = 1:2
                Matrix_Gex(x, y, z, j) = mean(Matrix_OGex(x, y, z, :, j));
                Matrix_G(x, y, z, j) = mean(Matrix_OG(x, y, z, :, j));
            end
        end
    end
end



%% Current Matrix Data Save
path = 'C:\Simulation_Data\NanoArray\';
Dn=dir(path);
Dnum=length(Dn);
for j= 1:5
    writematrix(Matrix_C(:,:,j), strcat(path,'#',num2str(Dnum),'Matrix_C_','.xls'),'sheet',j)
    writematrix(Matrix_Cex(:,:,j), strcat(path,'#',num2str(Dnum),'Matrix_Cex_','.xls'), 'sheet', j)
end

%% Current Matrix Graphics Output
for j = 1:5
    i = num2str(j);
    s = mesh(Matrix_C(:,:,j));
    s.FaceColor = 'flat';
    colorbar

    title(strcat('Current ligand#',i));
    ylabel(['Y dir xm1#', num2str(xm1)]);
    xlabel(['X dir ymax#', num2str(ymax)]);

    fig = gca;

    saveas(fig,strcat(path,'#',num2str(Dnum),'Current_',i,'.fig'));
end

for j = 1:5
    i = num2str(j);
    s = mesh(Matrix_Cex(:,:,j));
    s.FaceColor = 'flat';
    colorbar

    title(strcat('Current Ex ligand#',i));
    ylabel(['Y dir xm1#', num2str(xm1)]);
    xlabel(['X dir ymax#', num2str(ymax)]);

    fig = gca;

    saveas(fig,strcat(path,'#',num2str(Dnum),'Current_ex_',i,'.fig'));

end

pause(0.1)

%% Ligand Length Matrix Data Save
path = 'C:\Simulation_Data\NanoArray\';
Dn=dir(path);
Dnum=length(Dn);
for j= 1:5
    writematrix(Matrix_L(:,:,j), strcat(path,'#',num2str(Dnum),'Matrix_L_','.xls'),'sheet',j)
    writematrix(Matrix_Lex(:,:,j), strcat(path,'#',num2str(Dnum),'Matrix_Lex_','.xls'), 'sheet', j)
end

%% Ligand Length Matrix Graphics Output
for j = 1:5
    i = num2str(j);
    s = mesh(Matrix_L(:,:,j));
    s.FaceColor = 'flat';
    colorbar

    title(strcat('Ligand length ligand#',i));
    ylabel(['Y dir xml #', num2str(xm1)]);
    xlabel(['X dir ymax #', num2str(ymax)]);
    fig = gca;
    colormap(fig,jet);

    saveas(fig,strcat(path,'#',num2str(Dnum),'Ligand_',i,'.fig'));
end

for j = 1:5
    i = num2str(j);
    s = mesh(Matrix_Lex(:,:,j));
    s.FaceColor = 'flat';
    colorbar

    title(strcat('Ligand length Ex ligand#',i));
    ylabel(['Y dir xml#', num2str(xm1)]);
    xlabel(['X dir ymax#', num2str(ymax)]);
    fig = gca;
    colormap(fig,jet);

    saveas(fig,strcat(path,'#',num2str(Dnum),'Ligand_ex_',i,'.fig'));

end

pause(0.1)


%% Conductivity Matrix Data Save
path = 'C:\Simulation_Data\NanoArray\';
Dn=dir(path);
Dnum=length(Dn);
for z= 1:2
    for j= 1:5
        writematrix(Matrix_G(:,:,z,j), strcat(path,'#',num2str(Dnum),'Matrix_G_','with voltage_',num2str(z),'.xls'),'sheet',j)
        writematrix(Matrix_Gex(:,:,z,j), strcat(path,'#',num2str(Dnum),'Matrix_Gex_','with voltage_',num2str(z),'.xls'), 'sheet', j)
    end
end

%% Conductivity Matrix Graphics Output
for z = 1:2
    for j = 1:5
    i = num2str(j);
    s = mesh(Matrix_G(:,:,z,j));
    s.FaceColor = 'flat';
    colorbar

    title(strcat('Conductivity ligand#',i,' with voltage#',num2str(z)));
    ylabel(['Y dir xm1#', num2str(xm1)]);
    xlabel(['X dir ymax#', num2str(ymax)]);

    fig = gca;

    saveas(fig,strcat(path,'#',num2str(Dnum),'Conductivity_',i,'_with voltage_',num2str(z),'.fig'));
    end
end

for z = 1:2
    for j = 1:5
    i = num2str(j);
    s = mesh(Matrix_Gex(:,:,z,j));
    s.FaceColor = 'flat';
    colorbar

    title(strcat('Conductivity Ex ligand#',i,' with voltage#', num2str(z)));
    ylabel(['Y dir xm1#', num2str(xm1)]);
    xlabel(['X dir ymax#', num2str(ymax)]);

    fig = gca;

    saveas(fig,strcat(path,'#',num2str(Dnum),'Conductivity_Ex_',i,'_with voltage_',num2str(z),'.fig'));

    end
end
pause(0.1)

%% Data Save
output = struct('SmeGG', SmeGG,'HmeGG',HmeGG,'SmeGGStd', SmeGGStd,'HmeGGStd',HmeGGStd,'SmeG',SmeG,'HmeG',HmeG,'SmeGStd',SmeGStd,'HmeGStd',HmeGStd,'HNGG',HNGG,'HNG',HNG,'CNG',CNG,'CNGs',CNGs, 'Ea_alk', Ea_alk, 'Ea_ptc', Ea_ptc);
path = 'C:\Simulation_Data\NanoArray\';
filename = 'Randl';
Dn=dir(path);
Dnum=length(Dn);
pause(0.1)
savetime=strcat(date,'_',num2str(Dnum+1));
dlmwrite(strcat(path,filename,savetime,'.txt'), 'Vx','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), Vx,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'timelapse','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), timelapse,'-append','delimiter',',', 'precision',10,'newline','pc');

dlmwrite(strcat(path,filename,savetime,'.txt'), 'LM length','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), LM,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'LS Stdev','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), LS,'-append','delimiter',',', 'precision',10,'newline','pc');

dlmwrite(strcat(path,filename,savetime,'.txt'), 'NPgridrows','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), NPgridrows,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'Trials','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), Trials,'-append','delimiter',',', 'precision',10,'newline','pc');


dlmwrite(strcat(path,filename,savetime,'.txt'), 'SmeGG','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.SmeGG,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'HmeGG','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.HmeGG,'-append','delimiter',',', 'precision',10,'newline','pc');

dlmwrite(strcat(path,filename,savetime,'.txt'), 'SmeGGStd','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.SmeGGStd,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'HmeGGStd','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.HmeGGStd,'-append','delimiter',',', 'precision',10,'newline','pc');

dlmwrite(strcat(path,filename,savetime,'.txt'), 'SmeG','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.SmeG,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'HmeG','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.HmeG,'-append','delimiter',',', 'precision',10,'newline','pc');

dlmwrite(strcat(path,filename,savetime,'.txt'), 'SmeGStd','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.SmeGStd,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'HmeGStd','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.HmeGStd,'-append','delimiter',',', 'precision',10,'newline','pc');


dlmwrite(strcat(path,filename,savetime,'.txt'), 'HNGG','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.HNGG,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'HNG','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.HNG,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'CNG','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.CNG,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'CNGs','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.CNGs,'-append','delimiter',',', 'precision',10,'newline','pc');

dlmwrite(strcat(path,filename,savetime,'.txt'), 'Ea_alk','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.Ea_alk,'-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), 'Ea_ptc','-append','delimiter',',', 'precision',10,'newline','pc');
dlmwrite(strcat(path,filename,savetime,'.txt'), output.Ea_ptc,'-append','delimiter',',', 'precision',10,'newline','pc');


%% LigandEx Graphics Output
colorz=[0,0,0];
    for k=1:grids
        multiplier=0.99/grids;
        colorz(k,:)=[0,multiplier*k,multiplier*k];
    end
    figure;
    hold on

    for m=1:grids

        subplot(1,4,1)
        hold on
        semilogy(LM(1:5),HNGG(m,:),'-o', 'Color',colorz(m,:))

        subplot(1,4,2)
        hold on
        semilogy(LM(1:5),HNG(m,:),'-o', 'Color',colorz(m,:))

    end
    pause(0.1)

title(strcat('NProws', num2str(NPgridrows),', LMfactor=', num2str(lengthfactor),', LS=data, Vx=',num2str(Vx),' over ',num2str(Trials),' Trials Avging over ',num2str(grids),' Grids; toc=',num2str(timelapse)));
subplot(1,4,1)
xlabel 'Mean Length Angstroms'
ylabel 'Normalized Successful Hops'
legend('HNGG, darker = smaller Vx');
set (gca, 'yscale','log')

subplot(1,4,2)
xlabel 'Mean Length Angstroms'
ylabel 'Normalized Successful Hops'
legend('HNG, darker = smaller Vx');
set (gca, 'yscale','log')

subplot(1,4,3)
semilogy(LM,CNG)
xlabel('Mean Length Angstroms');
ylabel('CNG - PTCDI/Alkane');

subplot(1,4,4)
errorbar(LM, mean(CNG), std(CNG)/sqrt(length(CNG)));
xlabel('Mean Length Angstroms');
ylabel('CNG - PTCDI/Alkane');

saveas(gca, strcat(path,'#',num2str(Dnum),mfilename,'.fig'), 'fig');

pause(0.1)

end


toc
