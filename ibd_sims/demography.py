import msprime
import stdpopsim

###### European-like demography
N_current = 500000  # Current effective population size
N_bottleneck = 2100  # Size during bottleneck
T_bottleneck = 1400  # Time of bottleneck (generations)
growth_rate = 0.004  # Growth rate per generation

euro_bottleneck = msprime.Demography()

# Add population with exponential growth
euro_bottleneck.add_population(
    name="pop_0",
    initial_size=N_current,
    growth_rate=growth_rate
)

# Add bottleneck event
euro_bottleneck.add_population_parameters_change(
    time=T_bottleneck,
    initial_size=N_bottleneck,
    growth_rate=0
)


##### Exponential
expon = msprime.Demography()

expon.add_population(
    name="pop_0",
    initial_size=1_000_000,
    growth_rate=0.0307
)

expon.add_population_parameters_change(
    time=150,
    growth_rate=0
)


##### Recent bottleneck
himba = msprime.Demography()

himba.add_population(
    name="pop_0",
    initial_size=1000,
    growth_rate=0.06
)

himba.add_population_parameters_change(
    time=12,
    growth_rate=-0.042
)

himba.add_population_parameters_change(
    time=86,
    growth_rate=0
)



##### Constant Ne
constant_Ne = msprime.Demography()

constant_Ne.add_population(
    name="pop_0",
    initial_size=10000,
    growth_rate=0
)

##### Constnat Ne 100k
constant_Ne100k = msprime.Demography()

constant_Ne100k.add_population(
    name="pop_0",
    initial_size=100000,
    growth_rate=0
)


###### 2-Ppop OOA

def create_demography():
    # Initialize demography with initial populations
    demography = msprime.Demography()
    
    # Add populations (going backwards in time, we start with the most recent state)
    demography.add_population(name="pop_0", initial_size=501436.3, growth_rate=0.0195)
    demography.add_population(name="AFR", initial_size=432124.6, growth_rate=0.0166)
    demography.set_symmetric_migration_rate(["AFR", "pop_0"], 2.5e-05)
    
    # Set the most ancient population parameters
    demography.add_population_parameters_change(
        time=205, initial_size=9279.2, growth_rate=0.00307, population="pop_0"
    )
    demography.add_population_parameters_change(
        time=205, initial_size=14474.0, growth_rate=0, population="AFR"
    )


        # At 920 generations ago: Migration rate and EUR population changes
    demography.add_symmetric_migration_rate_change(time=920, rate=0.00015, populations=["pop_0", "AFR"])

    demography.add_population_parameters_change(
        time=920, initial_size=1861.0, growth_rate=0, population="pop_0"
    )

    demography.add_population_parameters_change(
        time=920, initial_size=14474.0, growth_rate=0, population="AFR"
    )

    
    # At 2040 generations ago: Mass migration and migration rate change
    demography.add_mass_migration(
        time=2040, source="pop_0", dest="AFR", proportion=1.0
    )

    demography.add_symmetric_migration_rate_change(time=2040, rate=0, populations=["pop_0", "AFR"])

    # At 5920 generations ago: AFR population size change
    demography.add_population_parameters_change(
        time=5920, initial_size=7310, population="AFR"
    )

    
    return demography

# Create and debug the demography
ooa2 = create_demography()

##### Ashkenazi
species = stdpopsim.get_species("HomSap")
model = species.get_demographic_model("AshkSub_7G19")
ashkenazi = model.model
ashkenazi.populations[6].name = "pop0"
ashkenazi.populations[6].extra_metadata["id"] = "pop0"



