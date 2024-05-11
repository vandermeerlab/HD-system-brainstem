
%D:\Github\hexamaze_old_replay_anal\behaviour

this_path = 'D:\Data\Hexamaze\prot_7\2018-01-25\r16\';
this_ppm_path = paths.poly_path;

%% Check for existence of ppm info
this_ppm_info = [paths.poly_path, thisSdata.name, '_ppmdata.m'];
if isfile(this_ppm_info)
    % Load the info
    ppm_data = load(this_ppm_info);
    ppm_data = ppm_data.ppm_data;
else
    % Get the video corresponding to this recording (if exists...) and ask
    % the user to draw relevant lines to compute the ppm

    this_vid_fn = exp_spd.vid_filenames(strcmp(exp_spd.raw_filenames, ...
        thisSdata.name));
    [vname] = get_vname(this_data_path, this_vid_fn{1});
    store = 1;
    [maze_ppm, pla_ppm] = get_ppm_from_user(vname, params, store_path, store);
    
    %% Create the video object
    disp('Loading video...')
    vidObj = VideoReader(vname);
    % read first frame
    first_f = read(vidObj, 1);
    % make it lighter
    first_f = imlocalbrighten(first_f);
    % Clear the video object (not sure if useful)
    clear('vidObj');
    fig_name = ['ppm determination figure for session ' thisSdata.name];
    this_fig = figure('Name', fig_name);
    % plot the maze image from the video
    imagesc(first_f);
    %% Draw the maze length line
    % try to have user draw a line
    disp('Please draw vertical line from end of top box to end of bottom box');
    roi = drawline('StripeColor', 'b');
    % compute the length of the line - in pixels
    length_maze_pix = pdist2(roi.Position(1,:), roi.Position(2,:));
    % compute pixel per meter ratio:
    ppm_maze = length_maze_pix / (params.maze_length);
    
    fprintf(['Maze: Line of %.1f pixels corresponds to %.3f m, so a ppm ' ...
        'of %.1f pix/m\n'], ...
        length_maze_pix, params.maze_length, ppm_maze)
    
    %% Draw the platform length line
    disp('Please draw vertical line from opposite edges of the central platform')
    
    roi = drawline('StripeColor', 'r');
    % compute the length of the line - in pixels
    length_pla_pix = pdist2(roi.Position(1,:), roi.Position(2,:));
    % compute pixel per meter ratio:
    ppm_pla = length_pla_pix / (params.platform_diam);
    fprintf(['Platform: Line of %.1f pixels corresponds to %.3f m, so a ppm ' ...
        'of %.1f pix/m\n'], ...
        length_pla_pix, params.platform_diam, ppm_pla);
    
    % Store these values somewhere so we don't have to ask again! e.g with the
    % masks

    ppm_data = {};
    ppm_data.ppm_pla = ppm_pla;
    ppm_data.ppm_maze = ppm_maze;
    ppm_data.ppm_fig = this_fig; % try to store the figure so we can reopen it to check later?
    
    if store
        save(this_ppm_info, 'ppm_data');
    end
end


function [vname] = get_vname(data_path, this_vid_fn)
    % Check that the video file exists and possibly renames it if the
    % extension is missing
    if exist([data_path this_vid_fn], 'file')
        % rename the file with '.avi' but make sure the avi file doesn't exist
        % yet
        if exist([data_path this_vid_fn '.avi'], 'file')
            disp(['Arg - both avi and original video files exist for file ' this_vid_fn]);
            keyboard
        else
            movefile([data_path this_vid_fn], [data_path this_vid_fn '.avi'])
            disp(['Renamed file ' this_vid_fn ' into ' this_vid_fn '.avi'])
        end    
    elseif ~exist([data_path this_vid_fn '.avi'], 'file')
        disp(['Warning! video file not found, please check path: ' data_path this_vid_fn])
        return
    end
    vname = [data_path this_vid_fn '.avi'];% the video file name
end