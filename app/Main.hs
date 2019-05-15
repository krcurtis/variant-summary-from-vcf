--------------------------------------------------------------------------------
-- Copyright 2019 Fred Hutchinson Cancer Research Center

{-# LANGUAGE OverloadedStrings #-}

module Main where


import VCF

import Conduit
import Statistics.Sample

import qualified Data.ByteString.Char8 as B

import qualified Data.Vector.Unboxed as V

import Text.Read (readMaybe)
import Data.Maybe (fromJust, isJust)

main :: IO ()
main = do
  B.putStrLn (B.intercalate "\t" ["snp_name", "chromosome", "position", "ref_allele", "alt_allele", "R2", "CAF", "missing", "total"])
  runConduitRes $ (main_pipeline n compute_snp_summary) where n = 64


extract_field :: [Bool] -> B.ByteString -> B.ByteString
extract_field keep_selector entry = head kept_sections
  where
    sections = B.split ':' entry
    kept_sections = [ section | (keep, section) <- zip keep_selector sections, keep]
    

read_maybe_double :: String -> Maybe Double
read_maybe_double a = readMaybe a



compute_snp_summary :: VariantRow -> B.ByteString
compute_snp_summary row = text
  where
    format_fields = B.split ':' (vformat row)
    dosage_field_keep = [ (x == "DS")    | x <- format_fields ]
    dosage_text = map (extract_field dosage_field_keep) (vpayload row)

    maybe_dosages = map (read_maybe_double . B.unpack) dosage_text
    dosages = [ fromJust d | d <- maybe_dosages, isJust d]

    vec = V.fromList dosages
    caf = (mean vec) / 2.0
    sample_variance = varianceUnbiased vec
    binomial_variance = 2*caf*(1-caf)
    rsquared = sample_variance / binomial_variance
    rsquared_text = B.pack . show $ rsquared
    caf_text = B.pack . show $ caf
    missing_text = B.pack . show $ (length maybe_dosages - length dosages)
    total_text = B.pack . show $ (length maybe_dosages)
    
    gecco_snp_name = (vchrom row) <> ":" <> (vpos row) <> "_" <> (vref row) <> "/" <> (valt row)
    text = B.intercalate "\t" [ gecco_snp_name, vchrom row, vpos row, vref row, valt row, rsquared_text, caf_text, missing_text, total_text]

