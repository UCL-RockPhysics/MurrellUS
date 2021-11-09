const t_ref = datetime2unix(DateTime(1904,01,01,0,0,0))

"""
	getrootinfo(fid::HDF5.File)

Read groups and dataset names for random file access
Inputs:
* fid: hdf5 file reference

return info = Dict(:start_time, :groups, :datasets)
"""
function getrootinfo(fid::HDF5.File)
        start_time = DateTime(read_attribute(fid,"Date")*read_attribute(fid,"Time"),"dd/mm/yyyyHH:MM:SS")
        groups = keys(fid)
        datasets = []
        for g = 1:length(groups)
                push!(datasets,keys(fid[groups[g]]))
        end
        info = Dict(
        :start_time => start_time,
        :groups => groups,
        :datasets => datasets
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

Return trace
"""
function getdata(fid::HDF5.File,info,chan::Int64,Ind::Int64)
        dataset = fid[info[:groups][chan]][info[:datasets][chan][Ind]]
        # Read attributes from dataset
        # Averaging = read_attribute(dataset,"Averaging")
        delay_us = read_attribute(dataset,"Delay (us)")
        # Fpulse_MHz = read_attribute(dataset,"Fpulse (MHz)")
        # gain_dB = read_attribute(dataset,"Gain (dB)")
        # PRF_us = read_attribute(dataset,"PRF (us)")
        # stacks = read_attribute(dataset,"Stacks")
        # tloop_s = datetime2unix(DateTime(maximum([info[:datasets][chan][Ind]]),"yyyymmddHHMMSS"))
        # t_s = t_conv(read_attribute(dataset,"Time (s)"))
        # pulse_V = read_attribute(dataset,"Vpulse (v)")
        # width_us = read_attribute(dataset,"Width (us)")

        # read waveform data
        trace = read(dataset)
        # f₀ = length(trace)/width_us*1e6
        return trace, delay_us
end
"""
    gettime_us(fid::HDF5.File, info, chan; indices = 1:length(fid[info[:groups][chan]]))
Get timestamp of ultrasonic data from .hdf5 file and shift to Julia timebase
Inputs:
* fid: Reference to open .hdf5 file
* info: a dictionary containing group and dataset information
* chan: index of the channel
Return t_us
"""
function gettime_us(fid::HDF5.File, info, chan; indices = 1:length(fid[info[:groups][chan]]))
        g = fid[info[:groups][chan]]
        t_us = zeros(length(g))
        for k in indices
            dataset = g[info[:datasets][chan][k]]
            t_us[k] =  t_conv(read_attribute(dataset,"Time (s)"))
        end
        return t_us
end

"""
    t_conv(t_in)
Convert Labview time to unix time
Inputs:
* t_in: LabView timestamp in seconds
Return t_conv
"""
function t_conv(t_in)
        t_conv = t_in+t_ref
end

"""
    trackarrivals(fid, gropuname, master, t0, f0, fcc, front, back)

Track change in arrival times by cross-correlation with a master waveform.
Inputs:
 * fid: HDF5 file reference
 * info: dictionary containing file structure information
 * chan: Indice of the channel
 * master: a master waveform
 * t0: the arrival time in the master waveform (in seconds)
 * f0: the sampling frequency (must be the same for all)
 * fcc: the resampling frequency
 * front: size of front window
 * back: size of back window

Return DT
"""
function trackarrivals(fid::HDF5.File, info, chan, master, t0, args... ; indices = 1:length(fid[info[:groups][chan]]))
    DT = zeros(length(indices))
    t = t0
    for k in indices
        print
        trace,~ = getdata(fid, info, chan, k)
        DT[k], _ = finddt(master, t0, trace, t, args...)
        t = t0 - DT[k]
    end
    return DT
end

"""
    trackarrivals_hpass(fid, gropname, master, t0, f0, fcc, front, back)

Track change in arrival times by cross-correlation with a master waveform.
Inputs:
 * fid: HDF5 file reference
 * info: dictionary containing file structure information
 * chan: Indice of the channel
 * master: a master waveform
 * t0: the arrival time in the master waveform (in seconds)
 * fpass: the high pass frequency
 * f0: the sampling frequency (must be the same for all)
 * fcc: the resampling frequency
 * front: size of front window
 * back: size of back window

Return DT
"""
function trackarrivals_hpass(fid::HDF5.File, info, chan, master, t0, filter, args... ; indices = 1:length(fid[info[:groups][chan]]))
    DT = zeros(length(indices))
    t = t0
    for k in indices
        print
        trace,~ = getdata(fid, info, chan, k)
        traceF = filtfilt(filter,trace)
        DT[k], _ = finddt(master, t0, traceF, t, args...)
        t = t0 - DT[k]
    end
    return DT
end

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

    # indices of picks
    i1 = max(1, trunc(Int, t1*f0))
    i2 = max(1, trunc(Int, t2*f0))

    (imin, imax) = minmax(i1,i2)

    # ake safe windows
    N = length(s1)
    back_safe  = min(back, imin-1)
    front_safe = min(front, N-imax)

    # trim signals
    offset = imin-back_safe-1
    win = tukey(front_safe+back_safe+1,0.6)


    s1trimmed = s1[imin-back_safe:imax+front_safe]
    s1trimmed[1:i1-imin] .= 0
    s1trimmed[i1+front_safe+1-offset:end] .= 0
    s1trimmed[i1-back_safe-offset:i1+front_safe-offset] .*= win

    s2trimmed = s2[imin-back_safe:imax+front_safe]
    s2trimmed[1:i2-imin] .= 0
    s2trimmed[i2+front_safe+1-offset:end] .= 0
    s2trimmed[i2-back_safe-offset:i2+front_safe-offset] .*= win

    # time basis
    time = (0:back_safe+front_safe+imax-imin)/f0
    time_rs = collect(time[1]:(1.0/fcc):time[end])

    # interpolate
    s1_rs = cubicinterp(time, s1trimmed, time_rs)
    s2_rs = cubicinterp(time, s2trimmed, time_rs)

    # ccr
    lags = (-(size(s1_rs,1)-1):size(s1_rs,1)-1)
    ccr = crosscor(s2_rs, s1_rs, lags)

    # find the timeshift
    (X, I) = findmax(ccr)
    DT = lags[I]/fcc

    return DT, X
end

function plot_traces(fid::HDF5.File, info, chan; indices = 1:length(fid[info[:groups][chan]]))
    figure;
    trace,delay = getdata(fid, info, chan, 1)
    t = delay .+ (0.01:0.01:length(trace)/100)
    for k in indices
        trace,delay = getdata(fid, info, chan, k)
        plot(t,trace./2048 .+k, c="k")
    end
    return gcf()
end
