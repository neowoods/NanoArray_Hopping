%% Part 1 of Run: to check if new code works as well as the old
LM = [15.3408,16.2173,20.1895,20.8322,24.5774];
LS = [5.4829,3.8156,6.2199,4.7574,4.6362];
Vx=30;
Trials=20
tic
[na,SmeGC,HmeGC,~] = RandlGMPGClean(14,298,LM(1),25.5,17,LS(1),10.1,1.47,.98,.017,Vx,Trials);
toc 
SmeGC
HmeGC
tic
[na,SmeC,HmeC,~] = RandlLoopClean(14,298,LM(1),LS(1),.017,Vx,Trials);
toc
SmeC
HmeC
tic
[na,SmeG,HmeG,~] = RandlGMPGCM(14,298,LM(1),25.5,17,LS(1),10.1,1.47,.98,.017,Vx,Trials);
toc 
SmeG
HmeG
tic
[na,Sme,Hme,~] = RandlLoopCM(14,298,LM(1),LS(1),.017,Vx,Trials);
toc
Sme
Hme
HNGC = HmeGC./SmeGC;
HNC = HmeC./SmeC;
HNG = HmeG./SmeG;
HN = Hme./Sme;
CNC = HNGC/HNC
CN = HNG/HN

title 'Normalized Average Successful Hops over 10^4 Trials Averaging over Five Grids for Experimental Data'
xlabel 'Mean Length'
ylabel 'Normalized Successful Hops'
legend('PTCDI Transport Normalized to Alkane Using Normal Distribution for Length','PTCDI Transport Normalized to Alkane Using Lorentz Distribution for Length')
toc