#!/bin/bash
cd -- "$(dirname "$BASH_SOURCE")"

timestamp() {
  date +"%T"
}
timestamp

rm ./test/output/*
rm ./use_cases/01_huggett/01_baseline/output/*
rm ./use_cases/01_huggett/02_discrete_types/output/*
rm ./use_cases/02_aiyagari/04_non_convex_capital_tax/output/*
rm ./use_cases/03_two_assets_illiquid/01_baseline/output/*
rm ./use_cases/04_two_assets_risky/01_idiosyncratic_risk/output/*
rm ./use_cases/08_life_cycle/01_one_asset_life_cycle/output/*

matlab -nodisplay -batch "run test/run_all_tests"
matlab -nodisplay -batch "run use_cases/01_huggett/01_baseline/main"
matlab -nodisplay -batch "run use_cases/01_huggett/02_discrete_types/main"
matlab -nodisplay -batch "run use_cases/02_aiyagari/04_non_convex_capital_tax/main"
matlab -nodisplay -batch "run use_cases/03_two_assets_illiquid/01_baseline/main"
matlab -nodisplay -batch "run use_cases/04_two_assets_risky/01_idiosyncratic_risk/main"
matlab -nodisplay -batch "run use_cases/08_life_cycle/01_one_asset_life_cycle/main"

timestamp
