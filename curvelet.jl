function _beta(t)
    #While the absolute value is not included in any of the equations, it seems
    #to be required to properly calculate values with t~ -1.0
    if(abs(t) <= 1)
        sqrt(abs(sum([-5/32 21/32 -35/32 35/32 1/2].*(t.^[7 5 3 1 0]))))
    elseif(1 <= t)
        1.0
    else
        0.0
    end
end

const NU_A = 0.15 #constant from UDCT paper
const NU_B = 0.15

_w1p(t) = _beta((pi-abs(t))/(pi*NU_A))
_w0p(t) = _w1p(2*t*(1+NU_A))
_w0(w::(Float64,Float64)) = _w0p(w[1]).*_w0p(w[2])

_w1(w::(Float64, Float64)) = ((1-_w0(w)^2)^0.5)*_w1p(w[1])*_w1p(w[2])

_v1p(t::Float64,N) = _beta((2/N-1-t)/(2*NU_B/N))*_beta((t+1)/(2NU_B/N))

_vp(t,N,l) = _v1p(t-(2*(l-1))/N,N)

function T(theta)
    if(-pi/4 <= theta && theta <= pi/4)
        tan(theta)
    elseif(pi/4 < theta && theta < pi/2)
        2-tan(pi/2-theta)
    elseif(-pi/2 < theta && theta < -pi/4)
        -2-tan(pi/2-theta)
    elseif(theta < -pi/2)
        -2
    else#  theta > pi/2
        2
    end
end


function _vl(theta,N,l)
    if(l<N)
        _vp(float(T(theta)), N, l)
    else
        _vp(T(pi/2-theta),N,2*N+1-l)
    end
end

function _ulp(w::(FloatingPoint, FloatingPoint), N, l)
    _vl(atan2(w[2],w[1]),N,l)*_w1(w)
end


function _ul(w::(FloatingPoint, FloatingPoint), N::Int, l::Int)
    if(l==0)
        _w0(w)
    else
        result = 0
        a = 0
        b = 0
        for a=[-2pi 0.0 2pi]
            for b=[-2pi 0.0 2pi]
                result += _ulp((w[1]+a,w[2]+b), N, l)
            end
        end
        result
    end
end

function grid(x,size)
    #println("converting ", x, " to ", -pi+2.0*pi/size*(x-1))
    (-pi+2.0*pi/size*(x-1))
end

global asdf = 0
#require("memoize.jl")
#@memoize 
function generate_the_stupid_window(size::Int, N::Int, l::Int)
    fn = (x,y)->_ul((x,y),N,l)

    result = zeros(size, size)

    #Calculate the window, allowing it to have its full image expressed by
    #wraping around the edge
    for x=1:size
        for y=1:size
            result[y,x] = fn(grid(x,size), grid(y,size))
            if(y==7 && x==26)
                println("y=",y,",x=",x," results in ",fn(grid(x,size), grid(y,size)))
                global asdf = (grid(x,size),grid(y,size))
                println("^^^ from (",grid(x,size),",", grid(y,size),")")
            end
        end
    end
    return result
end

function teh_window(size, N, l)
    if(l==0)
        generate_the_stupid_window(size, N, l)*2
    else
        new_stupid(size, N, l)*2*sqrt(2)
    end
end

#performs the downsample by a factor of two in both directions
function downsample(X)
    fftshift(fft(ifft(ifftshift(X))[1:2:end,1:2:end]))
end

function downsample2(X)
    Xx = fftshift(X)
    y  = X[1:end/2,1:end/2]
    y += X[1:end/2,end/2+1:end]
    y += X[end/2+1:end,end/2+1:end]
    y += X[end/2+1:end,1:end/2]
    fftshift(y)/4
end

function upsample(X)
    N = size(X,1)
    xx = zeros(2*N,2*N)+0im
    xx[1:2:end,1:2:end] = ifft(ifftshift(X))
    fftshift(fft(xx))
end

function upsample2(X)
    Xx = fftshift(X)
    [Xx Xx; Xx Xx]
end


C1 = 6
C2 = 2*C1+1

function windows(Wsize,N)
    win = Array(Any, 2*N+1)
    for i=0:N
        print('.')
        win[i+1] = teh_window(Wsize,N,i)
    end
    for i=N+1:2N
        win[i+1] = flipud(win[i-N+1])'
        win[i+1][:,[2:end, 1]] = win[i+1][:, [1:end]]
    end
    win
end

transform(x) = transform(x, windows(size(x,1), C1))

function transform(x, win)
    println("FFT")
    S = fftshift(fft(x))

    tf = Array(Any, C2)
    print("Transforming")
    for i=1:C2
        print(".")
        tf[i] = S.*win[i]
    end
    println()

    ds = Array(Any, C2)
    print("Downsampling")
    for i=1:C2
        ds[i] = downsample2(tf[i])
        print(".")
    end
    println()
    ds
end

function restore(ds, win)
    Wsize = size(ds[1],1)*2
    
    us = Array(Any, C2)
    print("upsampling")
    for i=1:C2
        us[i] = upsample2(ds[i])
        print(".")
    end
    println()

    prev = zeros(Wsize,Wsize)
    print("Inverting")
    for i=1:C2
        print(".")
        prev += us[i].*win[i]
    end
    println()

    println("IFFT")
    println("imaginary : ",norm(imag(ifft(ifftshift(prev)))))
    real(ifft(ifftshift(prev)))
end


function new_mod(a::FloatingPoint,b::FloatingPoint)
    while(signbit(a)==1)
        a += b
    end
    mod(a,b)
end

function func_1(w1::Float64,w2::Float64)
    w1=new_mod(w1,2pi)
    (w1 >= (pi/2)*(1-NU_A)/(1+NU_A)) & (w1 <= pi*(1+NU_A)) & (abs(w2) < pi/2*(1+2*NU_B))
end

function interesting_region(size,N,l)

    result = zeros(size, size)

    #grid = linspace(-pi,pi,size)
    grid = -pi+2.0*pi/size*(0:size-1)

    #Calculate the window, allowing it to have its full image expressed by
    #wraping around the edge
    for x=1:size
        for y=1:size
            result[y,x] = func_1(grid[x], grid[y])
        end
    end
    return result
end

function refit(f, low, high)
    while(f(low-1) != 0)
        low -= 1
    end
    while(f(low) == 0 && low != high)
        low += 1
    end
    while(f(high+1) != 0)
        high += 1
    end
    while(f(high) == 0 && high != low)
        high -= 1
    end

    return (low, high)
end

#extra parameters are used to avoid some WEIRD compiler bug with reguards to
#lambda precision...
function run_cols(f, x, b, A, dir, N, l)
    #println("working on good ol' col ", mod(x-1,size(A,1))+1)
    #V = x::Float64->-pi+2.0*pi/size(A,1)*(x-1)
    M = size(A,1)
    (low, high) = refit(y->f(grid(x,M),grid(y,M)), b[1], b[2])
    if(mod(x-1,size(A,1))+1 == 2)
        #println("at the danger column...")
        #println("f(x,1)=>",f(x,1))
        #println("f(x,0)=>",f(x,0))
    end
    if(low == high)
        x_ = mod(x-1, size(A,1))+1
        y_ = mod(low-1,size(A,2))+1
        A[y_,x_] = _ul((grid(x_,M),grid(y_,M)), N, l)
        return
    end
    for y=low:high
        x_::Int = mod(x-1, size(A,1))+1
        y_::Int = mod(y-1, size(A,2))+1
        A[y_, x_] = _ul((grid(x_,M),grid(y_,M)), N, l)
    end
    run_cols(f, x+dir, (low, high), A, dir, N, l)
end


function new_stupid(size::Int, N::Int, l::Int)
    if(l>N)
        result = flipud(new_stupid(size,N,N-l))'
        result[:,[2:end, 1]] = result[:, [1:end]]
        return result
    end
    fn = (x,y)->_ul((x,y),N,l)

    result = zeros(size, size)

    #grid = linspace(-pi,pi,size)
    #grid = -pi+2*pi/size*(0:size-1)

    #Calculate the window, allowing it to have its full image expressed by
    #wraping around the edge
    #for x=1:size
    #    for y=1:size
    #        result[y,x] = fn(grid[x], grid[y])
    #    end
    #end

    V = x->-pi+2*pi/size*(x-1)
    #f = (x,y)->fn(idxToVal(x), idxToVal(y))

    #-pi/2..3*pi/2 FIXME see eqn 28. T(theta)
    interp = x->(-1/2+1*x/N)
    slope  = (atan(pi*interp(l-1))+atan(pi*interp(l)))/2
    #println("slope=",slope)

    x::Int = floor(size*4/5) #TODO enhance this
    y::Int = round((size/(2*pi)*(V(x)*slope+pi))) + 1

    if(fn(V(x),V(y)) == 0)
        println("Warning: bad y value(",y,"), correcting...")
        #Bad guess, which should not happen if the above math is
        #right
        #oh well, do a quick search to find a valid point
        for ny = 1:size
            if(fn(V(x),V(ny)) != 0)
                y = ny
                break
            end
        end
    end


    @assert fn(V(x),V(y)) != 0
    #println("idx is (",x,",",y,")")
    #println("testing at x=",idxToVal(x), " y=", idxToVal(y))
    #println("0 != ", f(x,y))

    bounds = refit(y->fn(V(x),V(y)), y, y)
    #println("bounds = ", bounds)
    run_cols(fn, x, (y,y), result, -1, N, l)
    run_cols(fn, x+1, (y,y), result, +1, N, l)

    #destroy any bad values
    return result

    #Try to define the windows in a minimal number of steps utilzing a semi 
    #generic algorithm as defining the ROS cleanly is a PITA
    #For l!=0 it will be effectively a blurred image of a wedge with the inner
    #lowpass filter removed, the bluring and the fact that this is along an
    #angle make this annoying to deal with, however asside from the concave
    #region along corners this is actually a fairly easy shape to navigate
    #through
    #The algorithm is as follows:
    #- if(l=0), use a very simple easy window
    #- otherwise:
    #- based upon l/N the approximate angle can be found
    #- follow the exterior of the low pass to find the start of the window's ROS
    #- find a convient starting point north south east west depending upon the
    #  angle
    #- generate columns/rows of output over continious ROS
    #- dance and hope the overhead of this method is not too much
end



using Images
lena = imread("peppers.gif")
println("data is ", size(lena))
sample_data = reshape(lena[1,:,:],512,512)'*1.0
#
##sample_data = rand(1024,1024)
win = windows(size(sample_data,1), C1)
##sample_data -= mean(sample_data)
###sample_data = data
@time tmp = transform(sample_data, win)
##ttmp = transform(ifft(ifftshift(tmp[1])))
###ttmp[1] = fftshift(fft(restore(transform(real(ifft(ifftshift(ttmp[1]))))))).*teh_window(32,3,0)
##println("mean is = ", mean(ttmp[1]))
##tmp_one_old = tmp[1]
###tmp[1] = fftshift(fft(restore(ttmp)))
@time sample_hat  = restore(tmp, win)
println("Done")
println("Error was: ", norm(sample_data-sample_hat))
#end
#@iprofile report
#@iprofile clear
