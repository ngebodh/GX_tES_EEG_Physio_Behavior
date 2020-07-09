function data = fillBadChannelsFast(data,badchannels,chanlocs,method,channel_type)

% data = fillBadchannels(data,badchannels,method,locationfile)

%

%

% data          : samples x channels

% badchannels   : channel number x 1

% method        : either 'interp' or 'zeros'

% locationfile  : only used for interp method (channel number, theta, rho, channel name)

%

% jmad/jenma 2018

%

if nargin < 1

    data = randn(1000,64);

end

if nargin < 2

    badchannels = randi(64,3,1);

end

if nargin < 3

    method = 'interp';

end

if nargin < 4
      load('D:\GX Project\Analysis\ANTEEG32.mat');
%     load('..\data\location_file\BioSemi64.mat');

end

 

% channels_no = find(arrayfun(@(s) strcmpi(s.type,channel_type),chanlocs));
% channels_no = find(arrayfun(@(s) strcmpi(s.labels,channel_type),chanlocs));
 channels_no= 1:size(chanlocs,2);

%find the channel number for the chanlocs struct indicating the bad channels

if ~isempty(badchannels) && ~isempty(chanlocs)

    if isstruct(badchannels)

        for iChannel = 1:length(badchannels)

            badchannels_no(iChannel) = find(arrayfun(@(s) strcmpi(s.labels,badchannels(iChannel).labels),chanlocs(channels_no)));

        end

    else

        badchannels_no = badchannels;

    end

else

    badchannels_no = [];

end

 

Nchannels = size(data,2);

Nsamples = size(data,1);

 
 method = 'interp';
if ~isempty(data) && ~isempty(badchannels)

    if strcmpi(method,'interp')

        fprintf('Starting interpolation...')

       

        %% read the location file

        %calcuate the X and Y coordinates i.e. the 3D coordinated projected in to 2D space

        theta = [chanlocs(channels_no).theta];

        theta_rad = (theta+90)/360*2*pi; % rotate the coordinate system and convert to radian

        rho  = [chanlocs(channels_no).radius]; %get the lenght of the vector

        [X,Y] = pol2cart(theta_rad,rho); %convert from polar to cartesian coordinate system

       

        goodchannels_no = setdiff(1:Nchannels,badchannels_no)';

        

        fprintf('%d ',length(badchannels_no)), fprintf('channel(s) '), fprintf('%s ',chanlocs(channels_no(badchannels_no)).labels)

       

        %good channels coordinates

        Xgc = X(goodchannels_no)';

        Ygc = Y(goodchannels_no)';

        

        %bad channels coordinates

        Xbc = X(badchannels_no)';

        Ybc = Y(badchannels_no)';

       

        % start interpolating

        parfor iSample = 1:Nsamples

            data_sample = data(iSample,:);

           

            if ~isempty(goodchannels_no)

                %get the eeg sample values for the good channels

                Vgc = data_sample(goodchannels_no)';

               

                %create a scatter interpolation function

                F = scatteredInterpolant(Xgc,Ygc,Vgc);

               

                %use the function to interpolate missing values

                data_sample(badchannels_no) = F(Xbc,Ybc);

            else

                data_sample(badchannels_no) = 0;

            end

            data(iSample,:) = data_sample;

        end

    elseif strcmpi(method,'zeros')

        fprintf('Filling with zeros...')

        fprintf('channel(s) '), fprintf('%d ',badchannels)

        data(:,badchannels) = 0;

    end

    fprintf('...done\n')

end