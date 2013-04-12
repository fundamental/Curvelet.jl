require("Profile")
using IProfile

#@iprofile begin
function magical_function(arg)
    #something to do with the radius?
    Ind_rad  = round(log2(max(arg)))
    idx::Array{Int} = 8*2.^floor(0.0:0.5:((Ind_rad-2)/2))
end

function curveletTransform(f::Matrix{Float64})
    (n,m) = size(f)
    Wsize = [n,m]

    F=fftshift(fft(f))

    #TODO consider fourier mask

    scales  = magical_function([n,m])
    Nscales = length(scales)

    coeff = Array(Any, Nscales+1)

    #lets assume no wrapping for now
    curvelet = get_curvelet([n,m])
    x = windowize(F.*curvelet, Wsize)
    low = ifft(ifftshift(x))

    @printf("Working on %d scales\n", Nscales)
    for scale=1:Nscales
        Nangles=scales[scale]
        coeff[scale] = Array(Any, Nangles)
        @printf("\nWorking on %d angles", Nangles)
        for angle=1:Nangles
            print('*')
            curvelet = get_curvelet([n,m], scale, angle)
            x = windowize(F.*curvelet, Wsize)
            coeff[scale][angle] = ifft(ifftshift(x))
        end
    end

    coeff[Nscales+1] = low
    coeff
end

function inverseCurveletTransform(coeff, dims)
    Nradial = length(coeff)-1
    low     = coeff[end]
    F       = fftshift(fft(low)).*get_curvelet(dims)

    for r=1:Nradial
        Nangle = length(coeff[r])
        for angle = 1:Nangle
            curvelet = get_curvelet(dims, r, angle)
            F += fftshift(fft(coeff[r][angle])).*curvelet
        end
    end
    return real(ifft(ifftshift(F)))
end

get_origin(dims) = [floor(dims[1]/2)+1, floor(dims[2]/2)+1]

get_curvelet(dims) = get_curvelet(dims, 0.0, 0.0)
function get_curvelet(dims, scale::Number, angle_::Number)
    (N, M) = dims
    origin = get_origin([N, M])

    map = magical_function([N,M])
    Nscales = length(map)

    #Get the low pass coefficients
    #These are composed of a small window around DC
    if(angle_==0)
        tmp = zeros(N,M)
        tmp[origin[1]-1:origin[1]+1,origin[2]-1:origin[2]+1] = [0.928241517645832 1 0.928241517645832;1 1 1;0.928241517645832 1 0.928241517645832]
        return tmp
    end

    @assert scale <= Nscales
    @assert angle_ <= map[scale]

    direction = (angle_-1)*2*pi/map[scale]+pi/(2*map[scale])

    (X, Y) = meshgrid(1-origin[1]:N-origin[1],1-origin[2]:M-origin[2])
    (theta, rad) = cart2pol(X,Y)
    rad=(2^(-1.0*scale)).*rad
    theta=(2^(ceil(scale/2)+1))*angle(exp(im*(theta+direction.*ones(N,M))))/pi

    if(scale==Nscales)
        return finalRadialWindow(rad).*angleWindow(theta)
    else
        return radialWindow(rad).*angleWindow(theta)
    end
end

#simple 2d meshgrid replacement
function meshgrid(x::Array{Float64},y::Array{Float64})
    return (repmat(x', length(y), 1), repmat(y, 1, length(x)))
end

meshgrid(x::Range1{Float64}, y::Range1{Float64}) = meshgrid([x], [y])

#restricted version of matlab's cart2pol
cart2pol(x, y) = (atan2(y,x), sqrt(x.^2+y.^2))

function finalRadialWindow(r)
    r  = abs(r)
    N = size(r)
    (r .>= 5/6) + cos((5*ones(N)-6*r)*pi/2).*((r .>= 2/3) & (r .< 5/6))
end

function angleWindow(theta)
    #t = abs(theta)
    N = size(theta)
    result = zeros(N)
    for i=1:size(theta,1)
        for j=1:size(theta,2)
            t = abs(theta[i,j])
            if(t <= 1/3)
                result[i,j] = 1
            elseif(t <= 2/3)
                result[i,j] = cos((3*t-1)*pi/2)
            end
        end
    end
    result
    #FIXME unneeded allocation is happening here
    #(t .<= 1/3) + (t .<= 2/3).*(t.>1/3).*cos((3*t-ones(N))*pi/2)
end

function radialWindow(r)
    r  = abs(r)
    N  = size(r)
    result = zeros(N)
    for i=1:size(r,1)
        for j=1:size(r,2)
            R = r[i,j]
            if(R>=5/6 && R<=4/3)
                result[i,j] = 1
            elseif(R>=2/3 && R<=5/6)
                result[i,j] = cos((5-6*R)*pi/2)
            elseif(R>4/3 && R<5/3)
                result[i,j] = cos((3*R-4)*pi/2)
            end
        end
    end
    result


    #FIXME unneeded allocation is happening here
    #x  = (r .>= 5/6) & (r .<= 4/3)
    #x += cos((5*ones(N)-6*r)*pi/2).*((r .>= 2/3) & (r .<= 5/6))
    #x += cos((3*r-4*ones(N))*pi/2).*((r .> 4/3)  & (r .<= 5/3))
    #x
end

#This should perform the wrapping portion of the windowing to get an equivilent
#window that can utilize a smaller FFT
function windowize(image, Wsize)
    return image
    (X,Y)  = Wsize
    result = 0*im+zeros(X,Y)
    origin = get_origin(Wsize)
    yo::Int = origin[1]-floor(Y/2)
    xo::Int = origin[2]-floor(X/2)
    xs::Int = size(image,1)
    ys::Int = size(image,2)

    #no shift
    xstart::Int = xo-ceil(xo/X)*X
    ystart::Int = yo-ceil(yo/Y)*Y

    xnow::Int    = xstart
    ynow::Int    = ystart
    ystart_::Int = ystart

    #TODO figure out what is happening here and why
    while(xnow<X)
        while(ynow<Y)
            #window is added, boundaries are taken care of

            x1 = max(2-ynow,1):min(Y,Y+ys-ynow-Y+1)
            y1 = max(2-xnow,1):min(X,X+xs-xnow-X+1)
            x3 = max(ynow,1):min(ynow+Y-1,ys)
            y3 = max(xnow,1):min(xnow+X-1,xs)
            if(length([x1]) != 0)
                result[x1,y1] += image[x3,y3];
            end
            ynow=ynow+Y;
        end
        ynow=ystart
        while(ynow>ystart_)
            ynow=ynow-Y;
        end
        while(ynow+Y<ystart_)
            ynow=ynow+Y;
        end
        ystart=ynow;
        xnow=xnow+X;
    end
    imagesc(abs(result))
    println("updates___")
    println(abs(result))
    sleep(2)
    result
end

#curveletTransform(ones(16,16))
#curveletTransform(ones(16,16))
#curveletTransform(ones(32,32))
#curveletTransform(ones(128,128))
#ct_res = curveletTransform(ones(8,8))
#end
#@printf("\n\n\n The report\n---------\n")
#@iprofile report
#@iprofile clear

#Curvelab
function fdct_wrapping_window(x)
    wr = zeros(size(x,1), size(x,2))
    wl = zeros(size(x,1), size(x,2))
    x[abs(x) .< 2^-52] = 0 #TODO unneeded?
    wr[(x .> 0)& (x .< 1)] = exp(1-1./(1-exp(1-1./x[(x .> 0)& (x .< 1)])))
    wr[x .<= 0] = 1
    wl[(x .> 0)& (x .< 1)] = exp(1-1./(1-exp(1-1./(1-x[(x .> 0)& (x .< 1)]))))
    wl[x .>= 1] = 1

    normalization = sqrt(wl.^2 + wr.^2)
    wr ./= normalization
    wl ./= normalization
    return (wr, wl)
end

#@iprofile begin
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
    #else
    #    generate_the_stupid_window(size, N, l)*2*sqrt(2)
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
    #N = int(size(X,1)/2)
    #low  = 1:N/2
    #mid  = N/2+1:N*3/2
    #high = N*3/2+1:2*N
    #result = zeros(N, N)
    #result += X[mid,mid]
    ##alias terms
    #result += [flipud(X[high,mid]);flipud(X[low,mid])]
    #result += [fliplr(X[mid,low]) fliplr(X[mid,high])]
    #result += [fliplr(flipud(X[high,low])) fliplr(flipud(X[high,high]));fliplr(flipud(X[low,low])) fliplr(flipud(X[low,high]))];
    #result
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
    #low  = [1:N/2]
    #mid  = [N/2+1:N*3/2]
    #high = [N*3/2+1:2*N]

    #Y = zeros(2*N,2*N)+0im
    #Y[N/2+1:N*3/2,N/2+1:N*3/2] = X
    #Y[[low,high],mid] = 2*X
    #Y[mid,[low,high]] = 2*X
    #corner_row = [low,low,high,high]
    #corner_col = [low,high,low,high]
    #for i=1:length(corner_row)
    #    Y[corner_row[i], corner_col[i]] = X[floor(i/N)+1,mod(i,N)+1]
    #end
    ##Y[,] = X
    #Y
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


function transform(x, win)
    N = size(x,1)
    println(N)
    

    #print("Generating filters")
    #win = Array(Any, C2)
    #for i=0:C2-1
    #    print('.')
    #    win[i+1] = teh_window(N,C1,i)
    #end
    #println()

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
    #get the 256 inner values?
    for i=1:C2
        ds[i] = downsample2(tf[i])
        print(".")
        #N = size(tf[i],1)
        #ds[i] = tf[i][N/4+1:N*3/4,N/4+1:N*3/4]
    end
    println()
    ds
end

function restore(ds, win)
    N = size(ds[1],1)*2
    
    #print("Generating filters")
    #win = Array(Any, C2)
    #for i=0:C2-1
    #    print('.')
    #    win[i+1] = teh_window(N,C1,i)
    #end
    #println()

    us = Array(Any, C2)
    print("upsampling")
    for i=1:C2
        us[i] = upsample2(ds[i])
        print(".")
        #N = size(ds[i],1)*2
        #us[i] = zeros(N,N)+0im
        #us[i][N/4+1:N*3/4,N/4+1:N*3/4] = ds[i]
    end
    println()

    prev = zeros(N,N)
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
#using Winston

#ww(Wsize,N,l) = generate_the_stupid_window(Wsize,N,l)
#N=128
#println("leftover is = ",sum(interesting_region(N,1,1).*ww(N,3,2)-ww(N,3,2)))
#println("speedup is x",(N*N*9)/sum(interesting_region(N,1,1)))
#println("best speedup of 1 is x",(N*N*9)/sum((ww(N,8,1).!=0)))
#println("best speedup of 2 is x",(N*N*9)/sum((ww(N,3,2).!=0)))
#println("best speedup of 3 is x",(N*N*9)/sum((ww(N,3,3).!=0)))
#println("best speedup of 4 is x",(N*N*9)/sum((ww(N,3,4).!=0)))
#println("best speedup of 5 is x",(N*N*9)/sum((ww(N,3,5).!=0)))
#println("best speedup of 6 is x",(N*N*9)/sum((ww(N,3,6).!=0)))
#sleep(20)

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

global wasd = 0

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
        if(x_ == 26 && y_ == 7)
            println("consider this mystery")
            println((x_,y_),(grid(x_,M),grid(y_,M)))
            println("...")
            println("the real answer is..., ", _ul((grid(x_,M),grid(y_,M)), 3,1))
            println("the real answer is..., ", _ul((grid(x_,M),grid(y_,M)), N,l))
            println("N=",N," l=", l)
            global wasd = (grid(x_,M), grid(y_,M))
        end
        A[y_, x_] = _ul((grid(x_,M),grid(y_,M)), N, l)
        println("It is done, ",A[7, 26])
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


function fdct(x)
    #perform a real to complex transform
    #Let the finest level contain curvelets
    #Lets assume the number of scales is ceil(log2(min(M,N)) - 3)
    #Finally lets have 8 angles taken at the second finest level
    #   __________
    #   |\ |   | /|
    #   |_\|___|/_|
    #   |  |   |  |  approximately
    #   |__|___| _|
    #   | /|   |\ |
    #   |/_|___|_\|
end



using Images
lena = imread("circle.png")
println("data is ", size(lena))
sample_data = reshape(lena[1,:,:],256,256)'*1.0
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
