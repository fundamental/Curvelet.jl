include("../src/Curvelet.jl")
#using Curvelet
using PyPlot

PyPlot.close("all")

# Load sample data

input = imread("test/peppers.gif")
#input = imread("test/lena.jpg")
#input = input[:,:,1]
sample_data = input*1.0
N = size(sample_data,1)

# Display unmodified data
figure(1)
gray()
imshow(sample_data)
title("Unmodified Data")


println("Performing Conversion...")
plan        = Curvelet.planCurveletTransform(N, Curvelet.C1)
tmp         = Curvelet.curveletTransform(sample_data, plan)
sample_hat  = Curvelet.inverseCurveletTransform(tmp, plan)
println("Done")
println()
println("(Small) Error from transform was: ", norm(sample_data-sample_hat))
println()

# Display a subset of coefficients
combo = zeros(size(tmp.HP[1]))
N = length(tmp.HP)
for i=1:N
    combo += abs.(tmp.HP[i])
end
figure(10)
title("High Frequency Curves")
imshow(log.(combo))

figure(11)
title("A Single Curvelet Band")
imshow(log.(abs.(tmp.HP[5])))

#
##Create an approximation
#
function vectorize_curvelet(cc::Matrix{Complex{Float64}})
    vec(cc)
end

function vectorize_curvelet(cc::Curvelet.CurveletCoeff)
    data = vectorize_curvelet(cc.LP)
    for i = cc.HP
        data = [data,vec(i)]
    end
    data
end

function build_sparse(image, plan, sparsity)
    tmp = Curvelet.curveletTransform(sample_data, plan)
    data = vectorize_curvelet(tmp)
    ordered = reverse(sort(abs.(data)))
    curvelet_threshold(tmp, ordered[sparsity])
    Curvelet.inverseCurveletTransform(tmp, plan)
end
data = Float64[0.0]
for i=tmp.HP
    data = vcat(data, vec(i))
end
#
data = vec(data)
ordered = reverse(sort(abs.(data)))

#threshold = ordered[30000]
threshold = ordered[300]

function curvelet_threshold(cc::Matrix{Complex{Float64}}, threshold)
    #cc[abs(cc) .< threshold] = 0
end

function curvelet_threshold(cc::Curvelet.CurveletCoeff, threshold)
    curvelet_threshold(cc.LP, threshold)
    for i = 1:length(cc.HP)
        cc.HP[i][abs.(cc.HP[i]) .< threshold] = 0
    end
end

function reconstruction_step(observation, sparsity)
    tmp = Curvelet.curveletTransform(observation, plan)
    data = Float64[0.0]
    for i=tmp.HP
        data = vcat(data, vec(i))
    end
    #
    data = vec(data)
    ordered = reverse(sort(abs.(data)))
    threshold = ordered[sparsity]

    curvelet_threshold(tmp, threshold)
    sample_hat  = Curvelet.inverseCurveletTransform(tmp, plan)
end

begin# reconstruction_example()
    #attempt a reconstruction problem
    #grab around 30% of the points
    px_rate = 0.30

    N = size(observation,1)
    good_px = zeros(size(sample_data))
    for i=1:N, j=1:N
        good_px[i,j] = (rand() < px_rate) ? 1 : 0
    end

    observation = sample_data.*good_px
    sample_hat  = copy(observation)
    
    figure(100)
    title("Approximation Under Missing Data - Missing")
    gray()
    imshow(abs.(observation))

    for iters=1:40
        tr = 1000

        if(iters > 10)
            tr = 10000
        elseif(iters > 20)
            tr = 40000
        end

        smp_next = reconstruction_step(sample_hat, tr)
        println("Appx error  was: ", norm(sample_data-smp_next), " @ iter ",
        iters)
        sample_hat[good_px.==0] = smp_next[good_px.==0]
        if(iters in [1, 5, 10])
            figure(iters+100)
            title("Intermediate Denoising Stage")
            imshow(sample_hat)
        end
    end

    println("Denoising Missing Data Example")
    println()
    println("Obs  error  was: ", norm(sample_data-observation))
    println("Appx error  was: ", norm(sample_data-sample_hat))
    println("Orig energy was: ", norm(sample_data))

    figure(151)
    title("Approximation Under Missing Data - Reconstruction")
    gray()
    imshow(abs.(sample_hat))
    figure(152)
    title("Approximation Under Missing Data - Diff")
    gray()
    imshow(abs.(sample_hat).-sample_data)
end


tmp = Curvelet.curveletTransform(sample_data, plan)
K = round(Int,length(sample_data)*0.10)
threshold = ordered[K]
curvelet_threshold(tmp,threshold)


sample_hat  = Curvelet.inverseCurveletTransform(tmp, plan)
println()
println()
println("Creating a sparse represntation")
println()
println("Sparsity:       ", 100*K/length(sample_data), "%")
println("total energy:   ", norm(sample_data))
println("Appx error was: ", norm(sample_data-sample_hat))

figure(200)
title("Approximation Under Sparsity - Original")
gray()
imshow(abs.(sample_data))
figure(201)
title(string("Approximation Under Sparsity - ",K," sparse"))
gray()
imshow(abs.(sample_hat))
