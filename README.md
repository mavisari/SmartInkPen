# Longitudinal study of the writing gesture in elderly subjects using a sensorized pen 
This repository contains the code developed for the analysis of handwriting signals acquired through a Smart Ink Pen in a longitudinal study involving elderly subjects. The aim of the project is to investigate whether changes in handwriting can be used as an early indicator of cognitive decline.

The project was developed within the context of a biomedical engineering thesis titled “Studio longitudinale del gesto di scrittura in soggetti anziani tramite penna sensorizzata” at the Politecnico di Milano.

The study focuses on analyzing handwriting tasks performed by elderly individuals using a sensorized pen capable of recording motion and force signals during writing. These signals are processed to extract quantitative indicators describing handwriting performance across multiple domains such as time, force, frequency, and movement fluidity.

The extracted indicators are then used to:
- evaluate reliability of handwriting features
- compute the Minimal Detectable Change (MDC)
- detect potential subject-specific decline
- build machine learning models to classify cognitive status based on MMSE scores.

## Study Background 
With the progressive aging of the global population, the early detection of cognitive decline and neurodegenerative diseases has become increasingly important. Traditional clinical evaluations rely mainly on subjective scales such as the Mini Mental State Examination (MMSE). However, these approaches may not capture subtle changes occurring in the early stages of cognitive decline.

Handwriting represents a complex motor and cognitive task involving:
- fine motor control
- visuomotor coordination
- motor planning
- cognitive processing

Because of this complexity, handwriting alterations can reflect early neurological changes. The use of a Smart Ink Pen allows collecting objective quantitative measurements during natural writing tasks, enabling ecological monitoring of motor and cognitive performance.

## Dataset 
The dataset used in this project contains recordings collected from 57 elderly participants who used the Smart Ink Pen in their home environment over a 10-week monitoring period. Participants were asked to perform standardized tasks including: 
- Archimedean spiral drawing
- sentence copying in cursive
- sentence copying in block letters
In addition to handwriting data, MMSE scores were collected before and after the monitoring period to evaluate possible changes in cognitive status.
The pen records multiple signals such as:
- linear acceleration
- angular velocity
- magnetic field
- force applied to the pen tip

These signals are used to compute handwriting indicators in the following domains:
- time
- force
- fluidity
- inclination
- frequency

## MATLAB SCRITPS 
The MATLAB scripts implement the full pipeline for signal analysis and statistical evaluation:
- ```matlab FunctionIndicatorsComputationTA ``` computes handwriting indicators from raw sensor signals collected by the Smart Ink Pen.
- ```matlab ScriptAnalisiAffidabilita.m ``` performs statistical analysis to evaluate the test-retest reliability of the extracted indicators.

