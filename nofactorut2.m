
clear; clc;
load dates;
load names;
svf = load('svymeansVYTA2.txt');
svy = load('svymeansVYTA2.txt');
%svf(622,:) = [];
%svy(622,:) = [];
load ferrors6;

%type fbeta.csv
%a = readmatrix('fbeta2.csv');
%a(1,:) = [];
type ph1.csv
ph1 = readmatrix('ph1.csv');
ph1(1,:) = [];
type ph2.csv
ph2 = readmatrix('ph2.csv');
ph2(1,:) = [];
type ph3.csv
ph3 = readmatrix('ph3.csv');
ph3(1,:) = [];
type ph4.csv
ph4 = readmatrix('ph4.csv');
ph4(1,:) = [];
type ph5.csv
ph5 = readmatrix('ph5.csv');
ph5(1,:) = [];
type ph6.csv
ph6 = readmatrix('ph6.csv');
ph6(1,:) = [];
type ph7.csv
ph7 = readmatrix('ph7.csv');
ph7(1,:) = [];
type ph8.csv
ph8 = readmatrix('ph8.csv');
ph8(1,:) = [];
type ph9.csv
ph9 = readmatrix('ph9.csv');
ph9(1,:) = [];
type ph10.csv
ph10 = readmatrix('ph10.csv');
ph10(1,:) = [];
type ph11.csv
ph11 = readmatrix('ph11.csv');
ph11(1,:) = [];

ph = zeros(621,132,11);
ph(:,:,1) = ph1;
ph(:,:,2) = ph2;
ph(:,:,3) = ph3;
ph(:,:,4) = ph4;
ph(:,:,5) = ph5;
ph(:,:,6) = ph6;
ph(:,:,7) = ph7;
ph(:,:,8) = ph8;
ph(:,:,9) = ph9;
ph(:,:,10) = ph10;
ph(:,:,11) = ph11;



% Compute objects from predictors
h   = 12;
fb  = sparse(fbetas);
thf = [svf(1,:).*(1-svf(2,:));svf(2,:);svf(3,:).^2];
xf  = svf(4:end-3,:);
gf  = svf(end-3+1:end,:);
[evf] = compute_uf22(xf,thf,h);

evv = zeros(622,132,12);
for p = 1:12;
    evv(:,:,p) = cell2mat(evf(p));
end

ut2 = zeros(622,132,12);
ut2(:,:,1) = evv(:,:,1);
evv(622,:,:) = [];
ut2(622,:,:) = [];

for k = 2:12;
    for i = 1:621;
        for j = 1:132;
            ut2(i,j,k) = (ph(i,j,k-1)*(ut2(i,j,k).^0.5) + evv(i,j,k));
        end
    end
end


utcsa = squeeze(mean(sqrt(ut2),2));
utcsa = real(utcsa);











