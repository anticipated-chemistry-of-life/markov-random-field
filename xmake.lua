set_project("acol")
set_version("0.1")

add_rules("mode.debug", "mode.release")

set_languages("c++20")

set_policy("generator.compile_commands", true)
set_policy("build.release.strip", false) -- allows to see the function names in CPU profiling
add_rules("plugin.compile_commands.autoupdate", { outputdir = "build" })



-- =========================================================
-- Options
-- =========================================================

option("lotus", { default = true, defines = "USE_LOTUS", description = "Whether to use simple error model or lotus." })

add_requires("openmp", { system = false })
add_requires("zlib", "armadillo", "fmt", "fast_float", { system = false })
add_requires("gtest", { optional = true, system = false })
add_requires("nlohmann_json", { system = false })
add_requires("stduuid", { system = false })
add_requires("openssl", { system = false })


-- =========================================================
-- coretools
-- =========================================================

target("coretools")
on_config(function(target)
    import("devel.git")
    local git_version = git.lastcommit() or "N/A"
    git_version = git_version:gsub("%s+", "")
    target:add("defines", 'GITVERSION="' .. git_version .. '"')
end)
set_kind("static")
add_files("coretools/core/**.cpp")
add_headerfiles("coretools/core/**.h")

add_includedirs(
    "coretools/core/",
    { public = true }
)
add_packages(
    "fmt",
    "zlib",
    "armadillo",
    "fast_float",
    "openmp",
    { public = true, system = false }
)
add_cxxflags("-Wall", "-Wextra")

-- =========================================================
-- stattools
-- =========================================================

target("stattools")
set_kind("static")

add_files("stattools/core/**.cpp")
add_headerfiles("stattools/core/**.h")

add_includedirs(
    "stattools/core",
    { public = true }
)

add_deps("coretools")

add_cxxflags("-Wall", "-Wextra")

-- =========================================================
-- Main executable
-- =========================================================

target("acol")
set_kind("binary")
add_packages("openmp", "nlohmann_json", "stduuid", "openssl")
add_files("main.cpp")
add_files("src/**.cpp")
add_includedirs("src")
add_deps("coretools", "stattools")
add_frameworks("CoreFoundation", "Security")
add_defines("DEVTOOLS", "DEV_LOCATION")
add_options("lotus")
add_cxxflags("-Wall", "-Wextra", "-Werror", "-Wpedantic", "-Wuninitialized",
    "-Wreturn-local-addr",
    "-Warray-bounds",
    "-Wnull-dereference",
    "-Wdouble-promotion",
    "-Wformat=2",
    "-Wundef"
)
after_load(function(target)
    local pkg = target:pkg("openmp")
    if pkg then
        local linkdirs = pkg:get("linkdirs")
        if linkdirs then
            for _, linkdir in ipairs(linkdirs) do
                target:add("rpathdirs", linkdir)
            end
        end
    end
end)

-- =========================================================
-- Unit tests
-- =========================================================

target("acol_unitTests")
set_kind("binary")
set_default(false)

add_files("src/**.cpp")
add_files("tests/**.cpp")

add_includedirs("src", "tests")

add_packages("gtest_main", "openmp", "gtest")
add_links("gtest_main", "gtest")

add_deps("coretools", "stattools")

add_defines("CHECK_INTERVALS")

add_cxxflags("-Wall", "-Wextra")
