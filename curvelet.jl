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

function windowize(image, Wsize)
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
    result
end

#println("one")
#curveletTransform(ones(16,16))
#println("two")
#curveletTransform(ones(16,16))
#println("three")
#curveletTransform(ones(32,32))
#println("four")
#curveletTransform(ones(128,128))
ct_res = curveletTransform(ones(8,8))
#end
#@printf("\n\n\n The report\n---------\n")
#@iprofile report
#@iprofile clear


