using Markdown

@testset "Docstring rendering (Markdown corruption guard)" begin
    # Julia's Markdown.jl parser has a known, structure/length-dependent bug
    # (found repeatedly during the 2026-07 documentation audit): certain
    # constructs -- a backtick immediately followed by an apostrophe, a
    # literal `*` for multiplication inside a code span combined with dense
    # multi-underscore identifiers, or a bare (non-backticked) identifier
    # with several underscores -- can make the parser misread a later,
    # unrelated underscore as an emphasis delimiter. The visible symptom is
    # a source identifier's underscores silently rendering as asterisks
    # (`z_mid` -> `z*mid`), which is invisible unless the *rendered* output
    # is inspected -- editing the raw string and eyeballing it is not
    # sufficient, since the trigger depends on the surrounding text as a
    # whole, not just the edited span in isolation. This testset is the
    # permanent, automated form of the check that found and fixed four such
    # corruptions during Audit 2 (Batches 1-3a), one of which had shipped
    # silently for a full prior commit before being caught by chance.
    #
    # Detector A (hard-failing): for every binding with a real project
    # docstring, compare the count of every underscore-bearing identifier
    # token between the raw source text and the live `Base.doc()`-rendered
    # text. A real docstring never legitimately loses an identifier's
    # underscores on render, so any drop is corruption, full stop -- this
    # is a precise, zero-false-positive check (verified against every
    # binding in the module during the audit that established it).
    #
    # Detector B (soft, informational only): flags any bare asterisk glued
    # to a word/backtick character with no space, after stripping legitimate
    # `**bold**` markers. This catches the same corruption class from a
    # different angle, but has real false positives on legitimate content
    # (math notation like `O(n*m)`, single-word `*italic*`, `.* ` broadcast,
    # `bbox_*` wildcard mentions) -- it is not solid enough to hard-fail on,
    # but is retained as a defense-in-depth signal for corruption patterns
    # Detector A's identifier-count method might not happen to disturb.

    src_dir = joinpath(pkgdir(HomotopyGetsReal), "src")

    # Extract every `"""..."""` docstring from a file, keyed by the leading
    # identifier of its first non-blank content line (the signature line),
    # and separately record `@enum TypeName Member1 Member2 ...` /
    # `@enum TypeName begin ... end` member lists, since Julia's doc system
    # renders an enum member's `Base.doc()` as its *type's* docstring, not
    # a docstring of its own -- `Docs.meta` never has a direct entry for the
    # member, so without this, every enum member would be structurally
    # invisible to Detector A.
    function scan_file(path)
        lines = readlines(path)
        n = length(lines)
        file_docmap = Dict{String,String}()
        file_enum_members = Dict{String,String}()
        i = 1
        while i <= n
            if occursin("\"\"\"", lines[i]) && count("\"\"\"", lines[i]) == 1
                start = i
                j = i + 1
                while j <= n && !occursin("\"\"\"", lines[j])
                    j += 1
                end
                stop = j
                raw = join(lines[start+1:stop-1], "\n")
                sig = ""
                for k in (start+1):(stop-1)
                    s = strip(lines[k])
                    if !isempty(s)
                        sig = s
                        break
                    end
                end
                m = match(r"^[A-Za-z_][A-Za-z0-9_!]*", sig)
                if m !== nothing
                    name = m.match
                    file_docmap[name] = get(file_docmap, name, "") * "\n" * raw

                    k2 = stop + 1
                    while k2 <= n && isempty(strip(lines[k2]))
                        k2 += 1
                    end
                    if k2 <= n
                        mm = match(r"^@enum\s+([A-Za-z_][A-Za-z0-9_]*)\s+(.*)$", lines[k2])
                        if mm !== nothing && mm.captures[1] == name
                            rest = strip(mm.captures[2])
                            if rest == "begin"
                                k3 = k2 + 1
                                while k3 <= n && strip(lines[k3]) != "end"
                                    mem = strip(lines[k3])
                                    isempty(mem) || (file_enum_members[mem] = name)
                                    k3 += 1
                                end
                            else
                                for mem in split(rest)
                                    file_enum_members[mem] = name
                                end
                            end
                        end
                    end
                end
                i = stop + 1
                continue
            end
            i += 1
        end
        return file_docmap, file_enum_members
    end

    docmap = Dict{String,String}()
    enum_members = Dict{String,String}()
    for f in readdir(src_dir; join = true)
        endswith(f, ".jl") || continue
        dm, em = scan_file(f)
        merge!(docmap, dm)
        merge!(enum_members, em)
    end

    mod = HomotopyGetsReal
    names_to_check = Symbol[]
    for n in names(mod; all = true)
        n === :HomotopyGetsReal && continue
        startswith(string(n), "#") && continue
        isdefined(mod, n) || continue
        b = Docs.Binding(mod, n)
        (haskey(Docs.meta(mod), b) || Docs.doc(b) !== nothing) || continue
        push!(names_to_check, n)
    end

    function id_counts(text::AbstractString)
        counts = Dict{String,Int}()
        for m in eachmatch(r"\b[A-Za-z][A-Za-z0-9]*_[A-Za-z0-9_]*\b", text)
            counts[m.match] = get(counts, m.match, 0) + 1
        end
        return counts
    end

    function detector_a_failures(rendered::AbstractString, source::AbstractString)
        rc = id_counts(rendered)
        sc = id_counts(source)
        failures = String[]
        for (tok, cnt) in sc
            rcnt = get(rc, tok, 0)
            rcnt < cnt && push!(failures, "$tok ($cnt in source, $rcnt in rendered)")
        end
        return failures
    end

    function detector_b_hits(rendered::AbstractString)
        stripped = replace(rendered, r"\*\*[^\*]+\*\*" => "")
        [m.match for m in eachmatch(r"[A-Za-z0-9_.`]\*|\*[A-Za-z0-9_.`]", stripped)]
    end

    n_checked_direct = 0
    n_checked_enum = 0
    n_excluded = 0
    detector_b_report = Pair{String,Vector{SubString{String}}}[]

    @testset "Detector A: $n" for n in names_to_check
        key = string(n)
        b = Docs.Binding(mod, n)
        directly_documented = haskey(Docs.meta(mod), b)

        source_text = if directly_documented && haskey(docmap, key)
            n_checked_direct += 1
            docmap[key]
        elseif haskey(enum_members, key) && haskey(docmap, enum_members[key])
            n_checked_enum += 1
            docmap[enum_members[key]]
        else
            n_excluded += 1
            nothing
        end

        if source_text === nothing
            # No project docstring exists for this binding (verified: Base's
            # generic fallback docs render instead, e.g. Base.length or an
            # undocumented private helper) -- nothing of ours to corrupt.
            @test !directly_documented
        else
            rendered = string(Base.doc(getfield(mod, n)))
            failures = detector_a_failures(rendered, source_text)
            @test isempty(failures)
            if !isempty(failures)
                @error "Docstring rendering corruption detected" binding=n failures
            end
            hits = detector_b_hits(rendered)
            isempty(hits) || push!(detector_b_report, key => hits)
        end
    end

    # `Base.length(reg::VertexRegistry)` extends a Base function rather than
    # defining a new HomotopyGetsReal-owned binding, so it never appears in
    # `names(mod; all=true)` at all -- check its docstring directly via the
    # same Markdown round-trip, sourced straight from the raw text.
    if haskey(docmap, "Base")
        base_ext_source = docmap["Base"]
        parsed = Markdown.parse(base_ext_source)
        rendered = string(Markdown.plain(parsed))
        failures = detector_a_failures(rendered, base_ext_source)
        @testset "Detector A: Base.* extensions" begin
            @test isempty(failures)
        end
    end

    @info "Docstring rendering guard: $(n_checked_direct) directly-documented + $(n_checked_enum) enum-mapped bindings checked by Detector A; $(n_excluded) excluded (confirmed no project docstring)."
    if !isempty(detector_b_report)
        @info "Detector B (informational only, not a failure): possible stray asterisks -- verify manually if any binding here was recently edited" detector_b_report
    end
end
