load_model <- function(model_name="model_B"){
  ## ------------------------------------------------------------
  ## 1. Model Definitions (using rxode2)
  ## ------------------------------------------------------------
  require(rxode2)
  # -- Model 1: The original model, renamed to 'yang_model'
  if(model_name=="yang_model") return(rxode2({
    # Define model parameters with dummy initial values.
    # These will be overridden by the optimizer during fitting.
    kp    <- 0.1
    Kp    <- 100
    kdmax <- 0.1
    Kd50  <- 50
    nd    <- 4
    kbys  <- 0
    vglc  <- 0.01
    
    # Model equations
    Gp    = 0.5 * (G + sqrt(G^2 + 1e-9^2)) # ~max(G,0)
    mu    = kp * Gp / (Kp + Gp)
    kS    = kdmax * Kd50^nd / (Gp^nd + Kd50^nd)
    kB    = kbys * ND
    gate  = 1/(1 + exp(-G/1e-3)) # ~Heaviside(G)
    
    # Differential equations for the state variables
    d/dt(NL) = (mu - kS - kB) * NL
    d/dt(ND) = (kS + kB) * NL
    d/dt(G)  = -vglc * NL * gate
  })) 
  
  
  # -- Model 2: The new Hsieh et al. (2022) model 🔬
  if(model_name=="hsieh_model") return( rxode2({
    # Define model parameters
    Vmax <- 0.1  # Max substrate uptake rate
    Km   <- 100  # Michaelis-Menten constant for uptake
    m    <- 0.01 # Maintenance energy rate (mass substrate per mass biomass per time)
    Yxs  <- 0.5  # True growth yield (mass biomass per mass substrate)
    kmax <- 0.1  # Max death rate under starvation
    
    # Intermediate calculations
    Gp = 0.5 * (G + sqrt(G^2 + 1e-9^2)) # Differentiable version of max(G, 0)
    qs = Vmax * Gp / (Km + Gp)         # Specific substrate uptake rate
    
    # Core of Hsieh model: Net growth rate depends on maintenance surplus/deficit.
    # This uses a C-style ternary operator: (condition) ? (value_if_true) : (value_if_false)
    # rxode2 compiles this efficiently.
    mu_net = (qs >= m) * (Yxs * (qs - m)) + (qs < m) * (kmax * (qs - m) / m)
    
    # The death rate for ND accumulation is the positive part of the net death rate
    # k_death = -min(0, mu_net). This is 0 if mu_net > 0 (growth) and > 0 if mu_net < 0 (death).
    k_death = -min(0, mu_net)
    
    # Differential equations
    d/dt(NL) = mu_net * NL       # Live cells change based on the net growth rate
    d/dt(ND) = k_death * NL      # Dead cells accumulate when live cells die
    d/dt(G)  = -qs * NL         # Substrate is consumed for both growth and maintenance
  })) 
  
  
  if(model_name=="model_A") return(rxode2({
    ## parameters (defaults for fitting)
    kp    <- 4.227693e-02
    theta <- 4.630158e+03
    na    <- 7.429587e+00
    g50a  <- 6.670658e-01
    kd    <- 5.935262e-01
    nd    <- 6.114026e+00
    g50d  <- 1.237991e-03
    v     <- 1.799280e-05
    m     <- 2.321256e+00
    g50c  <- 2.270614e-02
    g50y  <- 3.341612e+03
    ky    <- 3.731199e-03
    
    ## growth term
    growth = kp * NL * (1 - NL/theta) * G^na / (G^na + g50a^na)
    ## death term (starvation + bystander)
    death  = kd * NL * (ky/kd * Y/(Y + g50y) + 1 - G^nd / (G^nd + g50d^nd))
    ## substrate consumption
    cons   = v  * NL * G^m  / (G^m  + g50c^m)
    
    ## ODEs
    d/dt(NL) = growth - death
    d/dt(ND) = death
    d/dt(G)  = -cons
    d/dt(Y)  = NL
  }))
  
  
  # -- Model B: constrained glucose consumption & confluency‑driven death --
  
  if(model_name=="model_B") return(rxode2({
    ## parameters (defaults for fitting)
    theta <- 1.586453e+04   # carrying capacity
    kp    <- 1.622025e-02   # proliferation rate
    kd    <- 7.569006e-02   # death rate
    kd2   <- 4.858519e-02   # confluency‑driven death
    g50a  <- 4.570815e-02   # half‑max for growth gating
    na    <- 5.912126e+00   # Hill coefficient for growth
    g50d  <- 3.475508e-03   # half‑max for death gating
    nd    <- 4.167656e+00   # Hill coefficient for death
    v1    <- 6.851104e-05   # consumption rate component A
    v2    <- 1.319974e-06   # consumption rate component D
    epsilon <- 1e-9 ## small constant enables 0 glucose simulation
    
    ## gating functions
    actA  = 1 / (1 + (g50a / (G+epsilon))^na)
    inhD  = 1 - 1 / (1 + (g50d / (G+epsilon))^nd)
    confl = kd2 * NL^2 / theta
    
    ## ODE system
    d/dt(NL) = kp * NL * (1 - NL/theta) * actA
    - kd * NL * inhD
    - confl
    
    d/dt(ND) = kd * NL * inhD
    + confl
    
    d/dt(G)  = - NL * (v1 * actA + v2 * (1 / (1 + (g50d / (G+epsilon))^nd))) / 2
  })) 
  
  ## intracellular glucose store model
  if(model_name=="intra_model") return(rxode2({
    ## -----------------------------------------------------------------------
    ## Model Parameters
    ## -----------------------------------------------------------------------
    # Glucose Uptake
    Vmax_uptake <- 0.0001    # mM/(#*hr)      ; max glucose uptake rate per cell
    Km_uptake   <- 0.1       # mM              ; half-saturation for uptake
    
    # Growth & Maintenance from Internal Store (R)
    kp          <- 1.6e-2    # 1/hr            ; max growth rate
    Y_xr        <- 500       # #/mM            ; yield of cells from internal resource R
    m_r         <- 1e-5      # mM/(#*hr)       ; maintenance cost per cell
    
    # Death Rates
    kdStarv     <- 8e-2      # 1/hr            ; max death rate from starvation
    kw          <- 5e-8      # 1/(mM*hr)       ; toxicity constant for waste
    
    # Waste Dynamics
    deltaW      <- 0         # 1/hr            ; waste removal/decay rate
    
    # Growth & Death Switches (based on internal store R)
    r_half_g    <- 0.001     # mM/#            ; per-cell store R for half-maximal growth
    nr_g        <- 4         # unitless        ; hill coef for growth switch
    r_half_d    <- 0.00002   # mM/#            ; per-cell store R for half-maximal death
    nr_d        <- 4         # unitless        ; hill coef for death switch
    
    ## -----------------------------------------------------------------------
    ## Intermediate Calculations & Differential Equations
    ## -----------------------------------------------------------------------
    # Per-cell internal store (average)
    r_cell <- max(0,R / (NL + 1e-9))
    
    # Key rates
    glucose_uptake_rate <- Vmax_uptake * (G / (Km_uptake + G)) 
    growth_rate <- kp * (r_cell^nr_g / (r_half_g^nr_g + r_cell^nr_g)) 
    starvation_death_rate <- kdStarv * (r_half_d^nr_d / (r_half_d^nr_d + r_cell^nr_d)) 
    waste_death_rate <- kw * W  # SIMPLIFIED: Linear toxicity
    
    # State variables
    d/dt(G) = -glucose_uptake_rate*NL
    
    d/dt(NL) = (growth_rate - starvation_death_rate - waste_death_rate)*NL
    
    d/dt(ND) = (starvation_death_rate + waste_death_rate)*NL
    
    d/dt(R) = NL*glucose_uptake_rate - NL*((1/Y_xr) * growth_rate + m_r + r_cell * (starvation_death_rate + waste_death_rate))
    
    d/dt(W) = (glucose_uptake_rate - deltaW * W)*NL # SIMPLIFIED: Waste integrates glucose use
    
  })) 
  stop(paste("no available model named",model_name))
}