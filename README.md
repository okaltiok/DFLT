# DFLT
Open access data and algorithms for RSS-based DFLT

1. Download the *.m files to a "dflt" folder on your computer, create a subfolder experiments and download the *.mat files and *.jpg file there.

2. Add the "dflt" folder and its subfolders to the Matlab path.

3. Open demo.m and run the code with the instructions commented in demo.m.

/////////////////////////////////////////////////////////////////////////

The flow of the algorithm is

    % initialize parameters

    % Perform N+1 EM iterations
    for i = 0:N

        % run filtering recursion (forward)
        for k = 1:K
            % compute one recursion of the tracking filter
            tracking_filter()

            % compute radio tomographic image
            rti()

            % call measurement selection unit
            measurement_selection()
        end

        % run smoothing recursion (backward)
        rtss()

        % compute parameter estimates using the EM algorithm
        em()
    end
