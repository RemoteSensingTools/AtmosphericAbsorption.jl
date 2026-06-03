#=
HITRAN API key management for the authenticated HITRANonline (HAPI2) endpoints — needed
only to fetch NON-Voigt parameters (Hartmann-Tran / speed-dependent / line-mixing). The
key is held in memory for the session and is NEVER written to disk or the repository;
users supply their own from their hitran.org profile.
=#

const _HITRAN_API_KEY = Ref{Union{Nothing,String}}(nothing)

"""
    activate_hitran!(key)

Store your HITRAN API key for this session (from your user profile at https://hitran.org).
Required for authenticated fetches of non-Voigt (Hartmann-Tran / speed-dependent / line-
mixing) parameters. The key lives in memory only — it is never written to disk or committed.
Alternatively, set the `HITRAN_API_KEY` environment variable.

```julia
activate_hitran!("xxxxxxxx-xxxx-xxxx-xxxx-xxxxxxxxxxxx")
```
"""
activate_hitran!(key::AbstractString) = (_HITRAN_API_KEY[] = String(key); nothing)

"""
    hitran_api_key() -> String

Resolve the HITRAN API key from `activate_hitran!` (preferred) or the `HITRAN_API_KEY`
environment variable. Throws a helpful error if neither is set.
"""
function hitran_api_key()
    k = _HITRAN_API_KEY[]
    k !== nothing && return k
    k = get(ENV, "HITRAN_API_KEY", nothing)
    (k !== nothing && !isempty(k)) && return k
    error("""
          No HITRAN API key found. Non-Voigt parameter fetches need one.
          Get your key from your profile at https://hitran.org, then either:
              activate_hitran!("<your-key>")
          or set the environment variable HITRAN_API_KEY before starting Julia.
          """)
end

"""Whether a HITRAN API key is available (without throwing)."""
has_hitran_api_key() =
    _HITRAN_API_KEY[] !== nothing || !isempty(get(ENV, "HITRAN_API_KEY", ""))
