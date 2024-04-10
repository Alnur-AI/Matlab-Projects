%% Subtasks 1-2
clc
clear
NIOT = readmatrix("data.xlsx", 'Range', 'E3:BF56');
TotalOutput = readmatrix("data.xlsx", 'Range', 'BO3:BO56');
FinalConsumption = readmatrix("data.xlsx", 'Range', 'BI3:BM56');
NIOT = NIOT.*1187.44;
TotalOutput = TotalOutput.*1187.44;
FinalConsumption = FinalConsumption.*1187.44;
%% Subtask 3
DirectCosts(1:22,1:22) = NIOT(1:22,1:22);
DirectCosts(1:22,23:53) = NIOT(1:22,24:54);
DirectCosts(23:53,1:22) = NIOT(24:54,1:22);
DirectCosts(23:53,23:53) = NIOT(24:54,24:54);
x(1:22) = TotalOutput(1:22);
x(23:53) = TotalOutput(24:54);
x = x';
A = DirectCosts./x;
Ax = A*x;
FinalConsumption = sum(FinalConsumption, 2);
w(1:22) = FinalConsumption(1:22);
w(23:53) = FinalConsumption(24:54);
w = w';
Axw = Ax+w;
flag = 1;
for i = 1:size(A,1)
   if(det(A(1:i,1:i)) == 0)
       flag = 0;
   end
end
if flag
    disp('Productive');
else
    disp('Not Productive');
end
%% Subtask 4
[EigenVector,EigenNumber] = eig(A);
DEigenNumber = diag(EigenNumber);
[FPNumber,FPIndex] = max(DEigenNumber);
FPVector = EigenVector(1:size(EigenVector,2,FPIndex))';
Lambda = FPNumber.*eye(size(A));
test = (A-Lambda)*FPVector;
%% Subtask 5
clc
MacroBranches = [3, 1, 18, 1, 2, 1, 3, 5, 1, 4, 3, 1, 5, 1, 1, 1, 1, 1];
Indexes = cumsum(MacroBranches);
Included = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
MB = MacroBranches.*Included;
MB(MB == 0) = [];
ID2 = Indexes.*Included;
ID2(ID2 == 0) = [];
ID1 = ID2-MB+1;
AggrA = zeros(size(MB));
Aggrx = zeros(size(MB,1));
Aggrw = zeros(size(MB,1));
for l = 1:size(MB,2)
    for k = 1:size(MB,2)
        AdditDC = DirectCosts(ID1(k):ID2(k), ID1(l):ID2(l));
        AdditSum = sum(sum((AdditDC)));
        Additx = x(ID1(l):ID2(l));
        AggrA(k,l) = AdditSum;
        Aggrx(l) = sum(Additx);
    end
    Aggrw(l) = sum(w((ID1(l):ID2(l))));
end
Aggrx = Aggrx';
Aggrw = Aggrw';
% sum(sum(AggrA))
% sum(sum(DirectCosts))
% sum(x)
% sum(Aggrx)
% sum(w)
% sum(Aggrw)
B = AggrA./Aggrx;
Bx = B*Aggrx;
Bxw = Bx+Aggrw;

flag = 1;
for i = 1:size(B,1)
   if(det(B(1:i,1:i)) == 0)
       flag = 0;
   end
end
if flag
    disp('Productive');
else
    disp('Not Productive');
end

%% Subtask 6
[AggrEigenVector,AggrEigenNumber] = eig(B);
AggrDEigenNumber = diag(AggrEigenNumber);
[AggrFPNumber,AggrFPIndex] = max(AggrDEigenNumber);
AggrFPVector = AggrEigenVector(1:size(AggrEigenVector,2,AggrFPIndex))';
AggrLambda = AggrFPNumber.*eye(size(B));
Aggrtest = (B-AggrLambda)*AggrFPVector;





