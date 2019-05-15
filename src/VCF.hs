--------------------------------------------------------------------------------
-- Copyright 2019 Fred Hutchinson Cancer Research Center
-- Module for reading VCFs as bytestrings

{-# LANGUAGE OverloadedStrings #-}


module VCF where

import Conduit

import System.IO (isEOF, hIsEOF, stdin)
import Data.Attoparsec.ByteString hiding (try)
import Data.Attoparsec.ByteString.Char8 (char8, double)
import GHC.Word
import qualified Data.ByteString.Char8 as B
import Data.Bits.Utils (c2w8)  -- use char8 instead
import Data.List (intersperse)
-- import Debug.Trace (trace)
import Control.Parallel.Strategies (parMap, rparWith, rdeepseq)
import Data.List.Split (chunksOf)


import Control.Monad (unless, replicateM)
import System.IO.Error
import Control.Exception 

import qualified Data.Map.Strict as M



--------------------------------------------------------------------------------


data VariantRow = VariantRowC { vchrom :: B.ByteString
                              , vpos :: B.ByteString
                              , vid :: B.ByteString
                              , vref :: B.ByteString
                              , valt :: B.ByteString
                              , vqual :: B.ByteString
                              , vfilter :: B.ByteString
                              , vinfo :: B.ByteString
                              , vformat :: B.ByteString
                              , vpayload :: [ B.ByteString ]
                              }


--------------------------------------------------------------------------------

bshow :: VariantRow -> B.ByteString
bshow row = (B.append header payload)
      where
        header = B.concat [ vchrom row
                       , tab
                       , vpos row
                       , tab
                       , vid row
                       , tab
                       , vref row
                       , tab
                       , valt row
                       , tab
                       , vqual row
                       , tab
                       , vfilter row
                       , tab
                       , vinfo row
                       , tab
                       , vformat row
                       , tab
                       ]
        payload = B.concat (intersperse tab (vpayload row))
        tab = B.pack "\t"

instance Show VariantRow where
    show row  = show $ bshow row

parse_variant_row :: Parser VariantRow
parse_variant_row = VariantRowC
    <$> (takeWhile1 (tab8 /=))
    <*> next_entry
    <*> next_entry
    <*> next_entry
    <*> next_entry
    <*> next_entry
    <*> next_entry
    <*> next_entry
    <*> next_entry
    <*> next_payload

next_payload :: Parser [B.ByteString]
next_payload = do
    values <- many1 next_entry
    return values

skiptab :: Parser GHC.Word.Word8
skiptab = char8 '\t'

tab8 = c2w8 '\t'

--first_entry :: Parser B.ByteString
--first_entry = takeWhile1 (tab8 /=)

--next_entry :: Parser B.ByteString
next_entry = do
    skiptab
    value <- takeWhile1 (tab8 /=)
    return value


--is_format_desc :: B.ByteString -> Bool
is_format_desc row = case (parse_results) of
    Done _ _ -> True
    _ -> False
  where parse_results = parse (string "##FORMAT") row 

is_sample_line row = case (parse_results) of
    Done _ _ -> True
    _ -> False
  where parse_results = parse (string "#CHROM") row 


--is_variant_row :: B.ByteString -> Bool
is_variant_row row = letter == (c2w8 '1') ||
                     letter == (c2w8 '2') ||
                     letter == (c2w8 '3') ||
                     letter == (c2w8 '4') ||
                     letter == (c2w8 '5') ||
                     letter == (c2w8 '6') ||
                     letter == (c2w8 '7') ||
                     letter == (c2w8 '8') ||
                     letter == (c2w8 '9') ||
                     letter == (c2w8 '0') ||
                     letter == (c2w8 'X') ||
                     letter == (c2w8 'Y') ||
                     letter == (c2w8 'M')
  where letter = (c2w8 . B.head) row




--------------------------------------------------------------------------------

-- conduit stuff to help parse file

--main_pipeline :: MonadIO m => (M.Map B.ByteString B.ByteString)  -> (M.Map B.ByteString B.ByteString) -> ConduitM () o m ()
--main_pipeline croc_map reach_map = my_get_line .| chunk_lines .| (transform_lines_with_r2_merging croc_map reach_map) .| stdoutC

main_pipeline :: MonadIO m => Int -> (VariantRow -> B.ByteString) -> ConduitM () o m ()
main_pipeline n converter = my_get_line .| (chunk_n_lines n) .| (transform_lines_with converter) .| stdoutC



my_get_line :: MonadIO m => ConduitM () B.ByteString m ()
my_get_line = gloop
  where
    gloop = do
      eof <- liftIO $ hIsEOF stdin
      case eof of
        True -> return ()
        False -> do
          line <- liftIO $ B.getLine
          yield line
          gloop


chunk_lines :: Monad m => ConduitM B.ByteString [B.ByteString] m ()
chunk_lines = cloop
  where
    cloop = do
      possible_lines <- replicateM 16 await
      let results = get_valid_values possible_lines
      case results of
          [] -> return ()
          _ -> do
            yield results
            cloop

chunk_n_lines :: Monad m => Int -> ConduitM B.ByteString [B.ByteString] m ()
chunk_n_lines n = cloop
  where
    cloop = do
      possible_lines <- replicateM n await
      let results = get_valid_values possible_lines
      case results of
          [] -> return ()
          _ -> do
            yield results
            cloop


get_valid_values :: [Maybe a] -> [a]
get_valid_values [] = []
get_valid_values (Nothing : xs) = get_valid_values(xs)
get_valid_values (Just x : xs) = x : get_valid_values(xs)




transform_lines_with_r2_merging :: Monad m => (M.Map B.ByteString B.ByteString)  -> (M.Map B.ByteString B.ByteString) -> ConduitM [B.ByteString] B.ByteString m ()
transform_lines_with_r2_merging croc_map reach_map = tloop
  where
    tloop = do
      mx <- await
      case mx of
        Nothing -> return ()
        Just lines -> do
          let line_groups = parMap (rparWith rdeepseq) (transform_vcf_line croc_map reach_map) lines
          let line_list = concat line_groups
          let with_endings = (intersperse "\n" line_list) ++ [ "\n"]
          yield (B.concat with_endings)
          tloop

transform_lines_with :: Monad m => (VariantRow -> B.ByteString) -> ConduitM [B.ByteString] B.ByteString m ()
transform_lines_with converter = tloop
  where
    tloop = do
      mx <- await
      case mx of
        Nothing -> return ()
        Just lines -> do
          let line_groups = parMap (rparWith rdeepseq) (transform_variant_row_with converter) lines
          let line_list = concat line_groups
          let with_endings = (intersperse "\n" line_list) ++ [ "\n"]
          yield (B.concat with_endings)
          tloop



new_info_lines = [ "##INFO=<ID=NONREACHR2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy (R-square)\">"
                 , "##INFO=<ID=REACHR2,Number=1,Type=Float,Description=\"Estimated Imputation Accuracy (R-square)\">" ]


transform_vcf_line :: (M.Map B.ByteString B.ByteString) -> (M.Map B.ByteString B.ByteString) -> B.ByteString -> [B.ByteString]
transform_vcf_line croc_map reach_map input
  | is_variant_row input = [variant_add_r2 croc_map reach_map input]
  | is_format_desc input = [input]
  | is_sample_line input = new_info_lines ++ [input]
  | otherwise = [input]

transform_variant_row_with :: (VariantRow -> B.ByteString) -> B.ByteString -> [B.ByteString]
transform_variant_row_with converter input
  | is_format_desc input = []
  | is_sample_line input = []
  | is_variant_row input = [ converter v ] 
  | otherwise = []
  where {
    Partial continuation = parse parse_variant_row input ;
    Done _ v = continuation (B.pack "") }

variant_add_r2 croc_map reach_map row = (bshow . (extend_info croc_map reach_map)) v
    where Partial continuation = parse parse_variant_row row
          Done _ v = continuation (B.pack "")
          


extend_info:: (M.Map B.ByteString B.ByteString) -> (M.Map B.ByteString B.ByteString) -> VariantRow -> VariantRow
extend_info croc_map reach_map v = VariantRowC (vchrom v) (vpos v) (vid v) (vref v) (valt v) (vqual v) (vfilter v) revised_info (vformat v) (vpayload v)
    where
        position = vpos v
        croc_r2 = case (M.lookup position croc_map) of
            Nothing -> "NA"
            Just x -> x
        reach_r2 = case (M.lookup position reach_map) of
            Nothing -> "NA"
            Just x -> x
        newinfo_croc = B.intercalate "=" [ "NONREACHR2", croc_r2]
        newinfo_reach = B.intercalate "=" [ "REACHR2", reach_r2]
        current_info = vinfo v
        sections = case current_info of
            ""  -> [newinfo_croc, ";", newinfo_reach]
            "." -> [newinfo_croc, ";", newinfo_reach]
            _ -> [current_info, ";", newinfo_croc, ";", newinfo_reach]
        revised_info = B.concat sections
        

transform_info :: B.ByteString -> B.ByteString
transform_info input_text= B.intercalate ";" keep_info_fields
    where
        sections = B.split ';' input_text
        fields = [ x | x:xs <- map (B.split '=') sections]
        keep_sections = [ (field == "AC") || (field == "AN") || (field == "SF") | field <- fields]
        keep_info_fields = [ section | (section, keep) <- zip sections keep_sections, keep]


transform_payload :: [Bool] -> B.ByteString -> B.ByteString
transform_payload drop_sections payload = B.intercalate ":" kept_sections
    where
        sections = B.split ':' payload
        kept_sections = [ section | (drop, section) <- zip drop_sections sections, not drop]
   

transform_sample_line :: B.ByteString -> B.ByteString
transform_sample_line column_header_row = B.intercalate "\t" updated_columns
    where
        chrom : pos : idf : ref : als : qual : filterf : info : format : samples = B.split '\t' column_header_row
        revised_samples = [ x | (x:xs) <- map (B.split '_') samples]
        updated_columns = chrom : pos : idf : ref : als : qual : filterf : info : format : revised_samples
