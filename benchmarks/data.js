window.BENCHMARK_DATA = {
  "lastUpdate": 1712678170131,
  "repoUrl": "https://github.com/JuliaSurv/NetSurvival.jl",
  "entries": {
    "Benchmark": [
      {
        "commit": {
          "author": {
            "email": "oskar.laverny@univ-amu.fr",
            "name": "Oskar Laverny",
            "username": "lrnv"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "b9406c5b9b60c44fda78b5ecf9a0384803b6c416",
          "message": "oops removed github pages (#14)",
          "timestamp": "2024-04-09T12:19:57+02:00",
          "tree_id": "9ea2ef0a8b23baf63efb8273a125196dba9f6555",
          "url": "https://github.com/JuliaSurv/NetSurvival.jl/commit/b9406c5b9b60c44fda78b5ecf9a0384803b6c416"
        },
        "date": 1712658342812,
        "tool": "julia",
        "benches": [
          {
            "name": "PoharPerme/colrec x slopop - explicit",
            "value": 229854922,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2339216\nallocs=29975\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "PoharPerme/colrec x frpop - formula",
            "value": 222140374,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6252856\nallocs=54870\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "GraffeoTest/colrec x slopop - explicit",
            "value": 228864155,
            "unit": "ns",
            "extra": "gctime=645968\nmemory=45582048\nallocs=31216\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      },
      {
        "commit": {
          "author": {
            "email": "oskar.laverny@univ-amu.fr",
            "name": "Oskar Laverny",
            "username": "lrnv"
          },
          "committer": {
            "email": "noreply@github.com",
            "name": "GitHub",
            "username": "web-flow"
          },
          "distinct": true,
          "id": "2c4e7f5f64acdf601493f6a7903d568b7df87a74",
          "message": "Stop checking on nightly temporarily (#16)",
          "timestamp": "2024-04-09T17:54:38+02:00",
          "tree_id": "af3cc9e30ed54f6e96171ac4e54c84329133e4c6",
          "url": "https://github.com/JuliaSurv/NetSurvival.jl/commit/2c4e7f5f64acdf601493f6a7903d568b7df87a74"
        },
        "date": 1712678169399,
        "tool": "julia",
        "benches": [
          {
            "name": "PoharPerme/colrec x slopop - explicit",
            "value": 219289399.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=2339216\nallocs=29975\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "PoharPerme/colrec x frpop - formula",
            "value": 220203662.5,
            "unit": "ns",
            "extra": "gctime=0\nmemory=6252856\nallocs=54870\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          },
          {
            "name": "GraffeoTest/colrec x slopop - explicit",
            "value": 229571200,
            "unit": "ns",
            "extra": "gctime=568013\nmemory=45582048\nallocs=31216\nparams={\"gctrial\":true,\"time_tolerance\":0.05,\"evals_set\":false,\"samples\":10000,\"evals\":1,\"gcsample\":false,\"seconds\":5,\"overhead\":0,\"memory_tolerance\":0.01}"
          }
        ]
      }
    ]
  }
}