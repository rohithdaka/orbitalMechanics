# a -> Semi-major Axis
# e -> eccentricity
# i -> inclination
# Ω -> Right Ascension of the Ascending Node
# ω -> Argument of Periapsis
# μ -> graviational parameter of the celestial body
include("Rot.jl")
include("NewtonRaphson.jl")

function CarVec(a,e,i,Ω,ω,M,μ)
    I = [1 0 0] #unit vectors 
    J = [0 1 0]
    K = [0 0 1]  
    (E ,I) = NewtonRaphson(M,e,M,1e-7,1) # Evaluating E  
    θ = 2*atan(sqrt((1+e)/(1-e))*tan(E/2))    
    
    if ( E < π && θ < 0) || (E > π && θ >0)  #Quadrant Check. 
        θ = θ + 2*π
    end  
    r = a*(1-e*cos(E)); #distance 
    v = sqrt(μ*(2/r - 1/a)) # velocity 
     
    R_polar = [ r 0 0] # in polar coordinates
    V_polar = [ sqrt(μ/a/(1-e^2))*e*sin(θ) sqrt(μ/a/(1-e^2))*(1+e*cos(θ)) 0] # vectors along polar directions 
     
    R = Rot(-Ω,3)*Rot(-i,1)*Rot(-θ-ω,3)*R_polar' #Transformations   
    V = Rot(-Ω,3)*Rot(-i,1)*Rot(-θ-ω,3)*V_polar'

    return (R,V)
end