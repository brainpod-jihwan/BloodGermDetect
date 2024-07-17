# BloodGermDetect
- https://img.shields.io/badge/any_text-you_like-blue

Blood Microbiome-Based Disease Diagnosis Pipeline

## 프로젝트 소개
최근 연구에 따르면, 폐혈증과 같은 미생물 감염병이 아니더라도 암과 같은 질환을 가진 환자의 혈액 내에서 존재해서는 안 되는 미생물의 DNA가 검출되고, 이것이 실제 질환과 연관이 있다는 보고가 등장했습니다. (Microbiome analyses of blood and tissues suggest cancer diagnostic approach; Poore et al. (Nature 2020)). 해당 연구에서는 혈액 내 저농도로 외래 미생물의 DNA가 검출되었으며, 이를 통해 환자의 질병을 예측하는 성과를 보였습니다.

특정 미생물의 존재 여부에 따라 환자의 질환을 예측할 수 있다는 것은, 조직 검사와 같은 침습적 검사 대신 획득이 쉬운 액체 생검만으로도 간단하게 질환을 진단할 수 있는 가능성을 열어줍니다.

본 프로젝트는 python과 기타 오픈소스를 활용하여 이러한 분석 파이프라인을 표준화하고, 발굴된 혈액 미생물에 대한 검출 검증을 위해 ddPCR 프라이머를 설계하여 제공합니다. 이를 기반으로 최종적으로 혈액 미생물 마커 기반 액상 생검 질환 진단의 기반 파이프라인을 제작하는 것을 목표로 합니다.

## 기능
- *WGS 데이터 분석*: 환자의 혈액 내 WGS 데이터를 인간 유전체에 매핑하고, 매핑되지 않은 리드(unmapped read)를 분석하여 인간 외래 미생물 DNA를 검출합니다.
- *ddPCR 프라이머 설계*: 검출된 미생물 DNA에 대한 ddPCR 프라이머를 설계하여 검출 민감도를 높입니다.
분석 파이프라인 표준화: 비효율적이고 어려운 기존 방법을 개선하여, 표준화된 분석 파이프라인을 제공합니다.

## 필요 조건
- Python 3.7 이상
- BioPython
- 기타 필요 패키지 (requirements.txt 참고)

## 참고 문헌
- Poore, G. D., et al. (2020). Microbiome analyses of blood and tissues suggest cancer diagnostic approach. Nature.
- Zmrzljak, U. P., et al. (2021). Detection of Somatic Mutations with ddPCR from Liquid Biopsy of Colorectal Cancer Patients. Genes.
- Ghezzi, H., et al. (2023). PUPpy: a primer design pipeline for substrain-level microbial detection and absolute quantification. doi:http://dx.doi.org/10.14288/1.0438913

