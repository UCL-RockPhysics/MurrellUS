# const t_ref = datetime2unix(DateTime(1904,01,01,0,0,0))

"""
	getrootinfo(fid::HDF5.File)

Read groups and dataset names for random file access, using dictionary database
Inputs:
* fid: hdf5 file reference

return info = Dict(:start_time, :groups, :datasets)
"""
function getrootinfo(fid::HDF5.File)
        start_time = DateTime(read_attribute(fid,"Date")*read_attribute(fid,"Time"),"dd/mm/yyyyHH:MM:SS") # read start time from file header
        groups = keys(fid) # get group keys (i.e. channel names) from .hdf5 file
        datasets = [] # create empty array to store dataset names
        for g = 1:length(groups) # loop through groups and extract dataset names using key function
                push!(datasets,keys(fid[groups[g]]))
        end
        info = Dict( # build dictionary to store data
        :start_time => start_time, # store start time
        :groups => groups, # store channel names
        :datasets => datasets # store dataset names
        )
end
"""
    getdata(fid::HDF5.File,info,chan::Int64,Ind::Int64)

Read a dataset and attributes from .hdf5 file using indices

Inputs:
* fid: .hdf5 file reference
* info: a dictionary containing group and dataset information
* chan: index of the channel
* Ind: the index of the dataset

Return trace, delay_us
"""
function getdata(fid::HDF5.File,info,chan::Int64,Ind::Int64)
        dataset = fid[info[:groups][chan]][info[:datasets][chan][Ind]] # open dataset from hdf5 file according to input channel and indice of dataset
        # Read attributes from dataset, most of these not required, but can be uncommented if other parameters need to be returned
        # Averaging = read_attribute(dataset,"Averaging")
        delay_us = read_attribute(dataset,"Delay (us)") # read delay attribute from individual data set
        # Fpulse_MHz = read_attribute(dataset,"Fpulse (MHz)")
        # gain_dB = read_attribute(dataset,"Gain (dB)")
        # PRF_us = read_attribute(dataset,"PRF (us)")
        # stacks = read_attribute(dataset,"Stacks")
        # tloop_s = datetime2unix(DateTime(maximum([info[:datasets][chan][Ind]]),"yyyymmddHHMMSS"))
        # t_s = t_conv(read_attribute(dataset,"Time (s)"))
        # pulse_V = read_attribute(dataset,"Vpulse (v)")
        # width_us = read_attribute(dataset,"Width (us)")

        # read waveform data
        trace = read(dataset) # read waveform data from dataset
        # f₀ = length(trace)/width_us*1e6
        return trace, delay_us # return the waveform trace and delay in microseconds
end
"""
    gettime_us(fid::HDF5.File, info, chan; indices = 1:length(fid[info[:groups][chan]]))
Get timestamp of ultrasonic data from .hdf5 file and shift to Julia timebase
Inputs:
* fid: Reference to open .hdf5 file
* info: a dictionary containing group and dataset information
* chan: index of the channel
Optional inputs:
* indices: an integer range corresponding to a subset of the data
Return t_us
"""
function gettime_us(fid::HDF5.File, info, chan; indices = 1:length(fid[info[:groups][chan]]))
        g = fid[info[:groups][chan]] # get name of group
        t_us = zeros(length(g)) # initialise array to be filled with survey timestamps
        for k in indices # loop through datasets
            dataset = g[info[:datasets][chan][k]] # read dataset into memory
            t_us[k] =  read_attribute(dataset,"Time (s)") # extract timestamp from open dataset
        end
        return t_us
end

# """
#     t_conv(t_in)
# Convert Labview time to unix time
# Inputs:
# * t_in: LabView timestamp in seconds
# Return t_conv
# """
# function t_conv(t_in)
#         t_conv = t_in+t_ref # convert LabView time to UTC by adding constant t_ref
# end

"""
    trackarrivals(fid, info, chan, master, t0, args...)

Track change in arrival times by cross-correlation with a master waveform.
Inputs:
 * fid: HDF5 file reference
 * info: dictionary containing file structure information
 * chan: Indice of the channel
 * master: a master waveform
Optional inputs
 * t0: the arrival time in the master waveform (in seconds)
 * f0: the sampling frequency (must be the same for all)
 * fcc: the resampling frequency
 * front: size of front window
 * back: size of back window
Keyword arguments
 * indices: an integer range corresponding to a subset of the data
Return DT
"""
function trackarrivals(fid::HDF5.File, info, chan, master, t0, args... ; indices = 1:length(fid[info[:groups][chan]]))
    DT = zeros(length(indices)) # get length of indice range of interest
    t = t0 # initialise variable t as the reference arrival time
    for k in indices # loop through indice range
        print # required to update variables(?)
        trace,~ = getdata(fid, info, chan, k) # get trace data, omit delay time
        DT[k-indices[1]+1], _ = finddt(master, t0, trace, t, args...) # find time shift between master and test waveform and store cumulative value relative to t0
        t = t0 - DT[k-indices[1]+1] # update time window center
    end
    return DT # return DT as a function output
end
#
# """
#     trackarrivals_hpass(fid, gropname, master, t0, f0, fcc, front, back)
#
# Track change in arrival times by cross-correlation with a master waveform.
# Inputs:
#  * fid: HDF5 file reference
#  * info: dictionary containing file structure information
#  * chan: Indice of the channel
#  * master: a master waveform
#  * t0: the arrival time in the master waveform (in seconds)
#  * fpass: the high pass frequency
#  * f0: the sampling frequency (must be the same for all)
#  * fcc: the resampling frequency
#  * front: size of front window
#  * back: size of back window
#
# Return DT
# """
# function trackarrivals_hpass(fid::HDF5.File, info, chan, master, t0, filter, args... ; indices = 1:length(fid[info[:groups][chan]]))
#     DT = zeros(length(indices))
#     t = t0
#     for k in indices
#         print
#         trace,~ = getdata(fid, info, chan, k)
#         traceF = filtfilt(filter,trace)
#         DT[k-indices[1]+1], _ = finddt(master, t0, traceF, t, args...)
#         t = t0 - DT[k-indices[1]+1]
#     end
#     return DT
# end

"""
    finddt(s1,t1,s2,t2,f0,fcc,front,back)
Compute time offset and maximum cross-correlation coefficient between two signals `s1` and `s2`, cut around times `t1` and `t2` with windoas of `front` points ahead and `back` points behind. The signals are resampled to frequency `fcc` (interpolated with cubic splines), multiplied with a tukey(x,0.6) window and padded with zeros. Need to have s1 and s2 of the same length. This is not checked so better be careful.
Input:
    s1, s2: input signals (same length !!)
    t1, t2: picks (in seconds)
    f0:     original sampling freq. (Hz)
    fcc:    resampling freq. (Hz)
    front:  nb of points in front window
    back:   nb of points in back window
Output:
    DT:     time difference for max ccr
    X:      value of max ccr
# Usage:
```
    (DT,X) = finddt(s1,t1,s2,t2,f0,fcc,front,back)
```
"""
function finddt(s1,t1,s2,t2,f0,fcc,front,back)
    # this was copied from AEprocessing, with no changes
    # indices of picks
    i1 = max(1, trunc(Int, t1*f0))
    i2 = max(1, trunc(Int, t2*f0))

    (imin, imax) = minmax(i1,i2) # to establish safe length of padding zeros

    # make safe windows
    N = length(s1) # get length of master waveform
    back_safe  = min(back, imin-1) # get minimum 'back' window i.e. after
    front_safe = min(front, N-imax) # get minimum 'front' window i.e. before

    # trim signals
    offset = imin-back_safe-1 # get offset betwen s1 and s2
    win = tukey(front_safe+back_safe+1,0.6) # compute tukey window


    s1trimmed = s1[imin-back_safe:imax+front_safe] # trim s1
    s1trimmed[1:i1-imin] .= 0 # pad front with zeros
    s1trimmed[i1+front_safe+1-offset:end] .= 0 # pad back with zeros
    s1trimmed[i1-back_safe-offset:i1+front_safe-offset] .*= win # window data

    s2trimmed = s2[imin-back_safe:imax+front_safe] # trim s2
    s2trimmed[1:i2-imin] .= 0 # pad front with zeros
    s2trimmed[i2+front_safe+1-offset:end] .= 0 # pad back with zeros
    s2trimmed[i2-back_safe-offset:i2+front_safe-offset] .*= win #window data

    # time basis
    time = (0:back_safe+front_safe+imax-imin)/f0 # create raw time vector
    time_rs = collect(time[1]:(1.0/fcc):time[end]) # create subsampled time vector

    # interpolate
    s1_rs = cubicinterp(time, s1trimmed, time_rs)
    s2_rs = cubicinterp(time, s2trimmed, time_rs)

    # ccr
    lags = (-(size(s1_rs,1)-1):size(s1_rs,1)-1) # create vector centered around zero of time lag
    ccr = crosscor(s2_rs, s1_rs, lags) # cross correlate the signals

    # find the timeshift
    (X, I) = findmax(ccr) # time shift corresponds to maxima of cross correlation function
    DT = lags[I]/fcc # compute timeshift in s

    return DT, X
end

"""
    plot_traces(fid::HDF5.File, info, chan; indices = 1:length(fid[info[:groups][chan]]), scale=0.1)
Plot traces in a 2d line chart to check arrival time computation and visualisation.
Input:
    fid:    reference to open .hdf5 file in memory
    info:   a dictionary containing information about the .hdf5 file
    chan:   the channel of data to plot
Keyword arguments:
    indices:indices or a range of traces to plot
    scale:  the normalisation scale of the data
Output:
    gcf():  handle of the created plot
# Usage:
```
    fig = plot_traces(fid, info, 3, indices = 1:10:1000, scale=0.1)
```
"""
function plot_traces(fid::HDF5.File, info, chan; indices = 1:length(fid[info[:groups][chan]]), scale=0.1)
    figure; # open new figure
    trace,delay = getdata(fid, info, chan, 1) # open the first waveform to create a time vector for plotting
    t = delay .+ (0.01:0.01:length(trace)/100) # create a time vector for plotting, in µs
    for k in indices
        trace,delay = getdata(fid, info, chan, k) # get data
        plot(t,scale.*trace.+k, lw=0.5, c="k") # plot data offsetting according to the trace number
    end
    return gcf()
end

"""
    V_calc!(P, info, range)
Compute normalised wave velocity history
Input:
    P: dictionary containing processed mechanical data
    info:   dictionary of information relating to mechanical data
    range: the portion of the experiment used to correct wavespeed with load
Output:
Adds the following parameters to the mechanical data dictionary
P[:DT]: the change in arrival time computed from cross correlation
P[:ΔV]: the normalised velocity change
"""
function V_calc!(P, info, range)
    ## Open waveform data
    pathUS = "/Users/christopherharbord/Dropbox/Research/UCL/raw_lab data/Murrell_US_data/"
    experiment=info[:name] # Get experiment name
    i = info[:I_US] # get indices of waveform data
    fid = h5open(pathUS*experiment*".h5") # Open .hdf5 file containing waveforms
    info_US = getrootinfo(fid) # Get information about file structure

    ## Mechanical data interpolation
    t_us = gettime_us(fid, info_US, 3) # get ultrasonic scan time for interpolation function
    P[:F_kN_i] = lininterp(P[:t_s],vec(P[:F_kN_j]), t_us) # interpolate corrected force
    P[:σ_MPa_i] = lininterp(P[:t_s],vec(P[:σ_MPa_j]), t_us) # interpolate differential stress
    P[:σ3_MPa_i] = lininterp(P[:t_s],vec(P[:Pc2_MPa]), t_us) # interpolate confining pressure
    P[:ε_i] = lininterp(P[:t_s],vec(P[:ε]), t_us) # interpolate sample axial strain

    ## Cross correlation algorithm
    master,delay = getdata(fid, info_US, 3, i[1]) # get the master waveform, corresponding to the experiment hit point
    P[:DT] = trackarrivals(fid, info_US, 3, master, info[:t1], 100e6, 1000e6, 50, 50, indices=i) # run the cross correlation algorithm to obtain arrival time

    ## Corrections for load and confining pressure
    P[:L_samp_m] = (1 .-(P[:ε_i][i] .-P[:ε_i][i[1]]))*info[:L_mm]*1e-3 # sample length during experiment
    # Compute length of pistons, subtracting confining pressure contribution and load contribution
    # this correction is computed using FEA and confirmed in experiment using a fused silica blank (run0146)
    P[:L_ass_m] = -(P[:σ3_MPa_i][i].*info[:U_P]) .- # Confining pressure contribution
                (P[:F_kN_i][i].*info[:U_F]) .+# Load contribution
                info[:L0] # initial length
    P[:T_ass_s] = P[:L_ass_m]./info[:V_ass] # Estimate travel time through pistons as a function of time
    P[:T_ini_s] = delay*1e-6 +info[:t0] .-P[:T_ass_s][1] # Estimate initial travel time through sample
    P[:V_ini] = P[:L_samp_m][1]/P[:T_ini_s] # estimate initial sample wavespeed

    ## Interface delay estimation
    x = P[:F_kN_i][i[range[1]]:i[range[end]]] # get force data to estimate interface delay
    y = -P[:DT][range] .-(P[:L_samp_m][range]./P[:V_ini] .+P[:T_ass_s][range])   # get delay time to estimate interface delay
    P[:m],_ = linfit(x,y) # get load dependant velocity change

    ## Final velocity calculations
    P[:T_samp] = delay*1e-6 +info[:t0] .-P[:T_ass_s] .-P[:DT] .-P[:F_kN_i][i]*P[:m] # correct travel time through sample for all delays
    V = P[:L_samp_m]./P[:T_samp] # compute velocity based on corrected travel time
    P[:V_ms] = V#./V[1] # store the velocity output
    P[:ΔV] = (V./V[1]).-1 # store the normalised velocity output
end
