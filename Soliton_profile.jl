function integrate(y, x) # trapez method
    so = y[1]*(x[2]-x[1]) + y[end]*(x[end]-x[end-1])
    for i in 2:length(y)-1
        so += y[i]*(x[i+1]-x[i-1])
    end
    0.5so
end;

# dr Needs to match resolution of soliton profile array file. Default = 0.00001
# Note that the spatial resolution of the profile must match 
# the specification of delta_x in the main code.

function soliton_profile_maker(dr = 1e-5)

    max_radius = 9.0
    n = round(Int64, max_radius/dr)
    g1(r, a, b, c) = 2*b*a - 2*c/r
    g2(r, a, d) = 4*π*a*a - 2*d/r
    
    s, beta, r90, diff = 0.0, 0.0, 0.0, 1.0
    optimised = false
    oo = 0
    la, lb, lc,  = zeros(n), zeros(n), zeros(n)
    ld, lr, intk = zeros(n), zeros(n), zeros(n)
    k1, k2 = zeros(4), zeros(4)
    k3, k4 = zeros(4), zeros(4);
    
    while optimised == false
        la[1] = 1.0
        lb[1] = s
        lc[1] = 0.0
        ld[1] = 0.0
        lr[1] = 0.001dr
        intk[1] = 0.0

        for i in 1:n-1

            k1[1] = lc[i]*dr
            k1[2] = ld[i]*dr
            k1[3] = g1(lr[i],la[i],lb[i],lc[i])*dr
            k1[4] = g2(lr[i],la[i],ld[i])*dr

            k2[1] = (lc[i] + 0.5k1[3])*dr
            k2[2] = (ld[i] + 0.5k1[4])*dr
            k2[3] = g1(lr[i]+0.5dr,la[i]+0.5k1[1],lb[i]+0.5k1[2],lc[i]+0.5k1[3])*dr
            k2[4] = g2(lr[i]+0.5dr,la[i]+0.5k1[1],ld[i]+0.5k1[4])*dr

            k3[1] = (lc[i] + 0.5k2[3])*dr
            k3[2] = (ld[i] + 0.5k2[4])*dr
            k3[3] = g1(lr[i]+0.5dr,la[i]+0.5k2[1],lb[i]+0.5k2[2],lc[i]+0.5k2[3])*dr
            k3[4] = g2(lr[i]+0.5dr,la[i]+0.5k2[1],ld[i]+0.5k2[4])*dr

            k4[1] = (lc[i] + k3[3])*dr
            k4[2] = (ld[i] + k3[4])*dr
            k4[3] = g1(lr[i]+dr,la[i]+k3[1],lb[i]+k3[2],lc[i]+k3[3])*dr
            k4[4] = g2(lr[i]+dr,la[i]+k3[1],ld[i]+k3[4])*dr

            la[i+1] = la[i]+(k1[1]+2k2[1]+2k3[1]+k4[1])/6
            lb[i+1] = lb[i]+(k1[2]+2k2[2]+2k3[2]+k4[2])/6
            lc[i+1] = lc[i]+(k1[3]+2k2[3]+2k3[3]+k4[3])/6
            ld[i+1] = ld[i]+(k1[4]+2k2[4]+2k3[4]+k4[4])/6
            lr[i+1] = lr[i]+dr
            intk[i+1] = (la[i]+(k1[1]+2k2[1]+2k3[1]+k4[1])/6)^2 * (lr[i]+dr)^2

            if (k1[1]+2k2[1]+2k3[1]+k4[1])/6.0 > 0.0
                s -= diff
                break
            end

            if la[i]+(k1[1]+2k2[1]+2k3[1]+k4[1])/6.0 < 0.0
                s += diff
                diff *= 0.1
                s -= diff
                break
            end

            if i == n-1
                optimised = true
                grad = (lb[i] - lb[i-1]) / dr
                co = lr[i]^2 * grad
                beta = lb[i] + co / lr[i]
            end
        end
        oo += 1
        if oo > 100
            break
        end
    end

    #Calculate full width at half maximum density:
    difflist = abs.(la.^2 .- 0.5)
    fwhm = 2lr[argmin(difflist)]
    #Calculate the (dimensionless) mass of the soliton:
    mass = integrate(intk,lr)*4pi
    #Calculate the radius containing 90% of the mass
    partial = 0.
    for i in 1:n
        partial = partial + intk[i]*4pi*dr
        if partial >= 0.9*mass
            r90 = lr[i]
            break
        end
    end

    partial = 0.
    for i in 1:n
        partial = partial + intk[i]*4pi*dr
        if lr[i] >= 0.5fwhm
            println("M_core = ", partial)
            break
        end
    end
    println("Full width at half maximum density is ", fwhm)
    println("Beta is ", beta)
    println("mass is ", mass)
    println("Radius at 90 mass is ", r90)
    
    lr, la, lb, beta, mass
end

function initsoliton!(psi_i, R, r_c, α, solprofile, δ_x)
    for index in CartesianIndices(psi_i)
        distfromcentre = hypot(
            R[index[1]] - r_c[1],
            R[index[2]] - r_c[2],
            R[index[3]] - r_c[3] )
        # Utilises soliton profile array out to dimensionless radius 5.6.
        if (sqrt(α) * distfromcentre <= 5.6)
            psi_i[index] = α * solprofile[round(Int,sqrt(α) * (distfromcentre / δ_x + 1))]
        else
            psi_i[index] = 0
        end
    end
end;