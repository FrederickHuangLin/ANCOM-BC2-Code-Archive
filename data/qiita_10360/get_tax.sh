{\rtf1\ansi\ansicpg1252\cocoartf2639
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 conda activate qiime2-2022.11\
\
qiime tools import \\\
  --input-path 131109_reference-hit.seqs.fa \\\
  --output-path qiita_10360_seqs.qza \\\
  --type 'FeatureData[Sequence]'\
\
qiime feature-classifier classify-sklearn \\\
  --i-classifier gg-13-8-99-515-806-nb-classifier.qza \\\
  --i-reads qiita_10360_seqs.qza \\\
  --o-classification qiita_10360_taxonomy.qza \\\
  --p-n-jobs 6 \
\
qiime tools export \\\
  --input-path qiita_10360_taxonomy.qza \\\
  --output-path taxonomy}