function good_data = mode_finder(harminv_dat, min_Q, min_amp, w_range, center, percent)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   mode_finder filters through harminv data and picks out possible good modes
%   given minimum Q values, amplitude values, and possibly a range of desired
%   frequencies
%
%   harminv_dat -> data in the format generated from "load_harminv_dat.m"
%   min_Q -> Minimum Q value (all modes with a Q lower then this will be filtered out)
%   min_amp -> Minimum amplitude (all modes with a lower amplitude will be filtered out)
%   w_range -> 1 or 0 value to specify if you want to also filter for a range of frequencies
%              if 1 it will take a center frequency and percentage (ratio from 0-1) and only 
%              return modes within center +- center*percent, if 0 the last to inputs dont matter
%   center -> Center frequency (in meep units) that you wish to filter modes for
%   percent -> A value from 0 to 1 (ratio) which will be multipled to the center
%              value to determine the range of frequencies which will be considered
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
good_ks = harminv_dat(((harminv_dat(:,4) > min_Q) & (harminv_dat(:,5) > min_amp)), 1);
good_freqs = harminv_dat(((harminv_dat(:,4) > min_Q) & (harminv_dat(:,5) > min_amp)), 2);
good_Qs = harminv_dat(((harminv_dat(:,4) > min_Q) & (harminv_dat(:,5) > min_amp)), 4);
good_amps = harminv_dat(((harminv_dat(:,4) > min_Q) & (harminv_dat(:,5) > min_amp)), 5);

if w_range % Are we specifiying a frequency range, if 0 the last two inputs dont matter
    good_ks = good_ks((good_freqs(:) > (center-center*percent)) & (good_freqs(:) < (center+center*percent)));
    good_freqs = good_freqs((good_freqs(:) > (center-center*percent)) & (good_freqs(:) < (center+center*percent)));
    good_Qs = good_Qs((good_freqs(:) > (center-center*percent)) & (good_freqs(:) < (center+center*percent)));
    good_amps = good_amps((good_freqs(:) > (center-center*percent)) & (good_freqs(:) < (center+center*percent)));
end

fprintf("k           freq            Q           amp\n")

for index = linspace(1,numel(good_ks),numel(good_ks))
    fprintf("%.3f %15.3e %*.0f %*.3e\n", good_ks(index), good_freqs(index), 6 + numel(num2str(round(good_Qs(index)))), good_Qs(index), 20 - numel(num2str(round(good_Qs(index)))),good_amps(index))
end

good_data = [good_ks good_freqs good_Qs good_amps];
end
