using AtmosphericAbsorption
using AtmosphericAbsorption.LineLists: resolve_molecule, resolve_isotopologue, molar_mass
using Test

const COPAR = joinpath(@__DIR__, "golden", "co_2100_2200.par")

@testset "species registry + resolvers" begin
    # molecule: symbol ≡ string ≡ integer; :ALL → -1; unknown throws
    @test molecule_number(:CO2) == 2 == molecule_number("CO2") == molecule_number(2)
    @test (molecule_symbol(2), molecule_symbol(5), molecule_symbol(30)) == (:CO2, :CO, :SF6)
    @test resolve_molecule(:ALL) == -1 && resolve_molecule(:CO) == 5 && resolve_molecule(7) == 7
    @test_throws ArgumentError molecule_number(:NotAMolecule)
    # isotopologue: :ALL/-1/:main/id and formula display
    @test resolve_isotopologue(2, :main) == 1 && resolve_isotopologue(2, :ALL) == -1
    @test resolve_isotopologue(:CO2, :main) == 1 && resolve_isotopologue("CO2", 3) == 3  # symbol/string mol
    @test resolve_isotopologue(2, 3) == 3
    @test isotopologue(2, 1) == "(12C)(16O)2"
    @test_throws ArgumentError resolve_isotopologue(2, :nope)
    # overridable mappings (use a free id so we don't clobber a canonical name)
    register_molecule!(:MyMol, 99)
    @test molecule_number(:MyMol) == 99 && molecule_symbol(99) == :MyMol
    register_isotopologue!(:CO2, Symbol("636"), 2)
    @test resolve_isotopologue(2, Symbol("636")) == 2
    # a name-only molecule (SF6 — .xsc-only) resolves + names, but has no line-list metadata
    @test_throws ArgumentError molar_mass(30, 1)
end

@testset "generic notation through the HITRAN port + attached partition" begin
    port = HitranPort(COPAR)
    db_sym, db_int = load_lines(port; mol = :CO), load_lines(port; mol = 5)
    @test 0 < length(db_sym) == length(db_int)          # :CO ≡ 5
    @test molecules(db_sym) == [:CO]
    @test db_sym.partition isa TIPS2021PF               # partition rides on the data
    @test LineByLineModel(db_sym; profile = Voigt()).partition isa TIPS2021PF  # model default
end

@testset "subset by molecule / mask preserves partition + meta" begin
    db = load_lines(HitranPort(COPAR); mol = :CO)
    sub = db[:CO]
    @test length(sub) == length(db) && sub.partition === db.partition && sub.meta === db.meta
    half = db[db.ν0 .< 2150]
    @test typeof(half) == typeof(db) && all(<(2150), half.ν0)
    @test eltype(db) == Float64
end

@testset "show: LineDatabase + LineByLineModel" begin
    db = load_lines(HitranPort(COPAR); mol = :CO)
    s = sprint(show, MIME"text/plain"(), db)
    @test all(n -> occursin(n, s), ("LineDatabase", "transitions", "TIPS-2021", "(12C)(16O)", "HITRAN"))
    m = LineByLineModel(db; profile = Voigt())
    sm = sprint(show, MIME"text/plain"(), m)
    @test all(n -> occursin(n, sm), ("LineByLineModel", "Voigt", "TIPS-2021", "CO", "HumlicekWeideman32"))
    @test occursin("$(length(db)) lines", sprint(show, db))       # compact form
end
