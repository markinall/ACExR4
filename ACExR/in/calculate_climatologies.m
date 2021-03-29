% calculate_climatologies

clear all

% specify input file
filename = 'ACmodel_rainfall_database.csv';

% sepcify years to derive climatology
% Use Year1=1970 and Year2 = 2006 to get overall climatology
% otherwise average over decades
Year1 = 1970;
Year2 = 2006;

% read data
data=dlmread(filename,',',3,0);
yy = data(:,1);
mm = data(:,2);
dd = data(:,3);
rain = data(:,5:end);
ncol = size(data,2);

% replace values of -999 with NaN
rain(rain==-999) = NaN;

% sepcify number of days per month
ndays_per_month = [31 28 31 30 31 30 31 31 30 31 30 31];

% Initialise output array
%ACrain = zeros(365,ncol);
ACrain = NaN(365,ncol);

% Select only the data for the required period
id = yy >= Year1 & yy <= Year2;
yy = yy(id);
mm = mm(id);
dd = dd(id);
rain = rain(id,:);

% derive climatology
yearday = 0;
for imn = 1:12
    for idy = 1:ndays_per_month(imn)
        yearday = yearday + 1;
        id1 = mm == imn;
        id2 = dd == idy;
        id = id1 & id2;
        ACrain(yearday,1) = Year1;
        ACrain(yearday,2) = imn;
        ACrain(yearday,3) = idy;
        ACrain(yearday,4) = yearday;
        for icol = 5:ncol
            %nrec = sum(rain(id,icol-4) > 0);
            rain2 = rain(id,icol-4);
            idreal = ~isnan(rain2);
            nrec = sum(idreal);
            if nrec > 0
                ACrain(yearday,icol) = sum(rain2(idreal)) / nrec;
            end
        end
    end
end
% Save data to file
dataout = [0 0 0 0 56.45 56.521 56.441 60.259 60.139];
dataout(2,:) = [0 0 0 0 -5.442 -5.323 -5.218 -1.287 -1.183];
dataout = [dataout; ACrain];
%fileout = ['Climatology_rainfall_',num2str(Year1),'-',num2str(Year2),'.csv'];
fileout = ['Test_rainfall_',num2str(Year1),'-',num2str(Year2),'.csv'];
dlmwrite(fileout,dataout,'delimiter',',','precision','%7.1f');
        