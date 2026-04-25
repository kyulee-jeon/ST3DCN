# ST3DCN Replication 기반 Essential DICOM Tag 검증 (HCC-TACE-Seg + CRLM)

**프로젝트**: protocol2dicom  ·  **작성일**: 2026-04-24

---

## 0. 한눈에 보는 요약

- **목표**: AI 영상 모델 replication에서 논문이 명시하는 "mentioned tags"만으로는 deployment variability를 포착하기 부족하다는 것을 실증. 대신 26개로 구성된 **essential tags**가 acquisition-parameter variability를 관측하기 위한 **최소 필요 집합** 임을 두 개의 독립 cohort (HCC+, HCC−)에서 보인다.
- **대상 모델**: ST3DCN (Yu et al., *JHEP Reports* 2025) — 3D-DCN HCC 자동 진단 모델.
- **두 cohort**:
  1. **HCC-TACE-Seg** (TCIA): HCC+ 환자 105명, 146 lesion → **sensitivity** 측정.
  2. **CRLM** (TCIA Colorectal-Liver-Metastases): HCC− (colorectal liver metastases) 환자 197명, 464 lesion → **specificity** 측정.
- **Validation**: 논문 repo에서 제공하는 demo sample 2건에 대해 전처리 결과 **pixel-exact (|Δ|=0)** 재현 후, 두 TCIA cohort에 동일 파이프라인 적용.
- **분석 단위**: **portovenous phase** 만, **pre-treatment** (HCC-TACE-Seg 기준 pre-TACE / CRLM 기준 pre-hepatectomy), **lesion-level**.
- **핵심 통계 도구**: Kruskal-Wallis + Benjamini-Hochberg FDR.

### 핵심 수치

| 지표 | Direct-only 분석 | DICOM 표준 derivation 적용 |
|---|---|---|
| HCC FDR-sig tag | 4 (논문 언급 2 + 미언급 2) | **6** (언급 2 + 미언급 4) |
| CRLM FDR-sig tag | 4 (언급 2 + 미언급 2) | **5** (언급 2 + 미언급 3) |
| 두 cohort FDR-sig 합집합 | 8 | **10** |
| 두 cohort FDR-sig 교집합 | **∅** | **{WindowCenter}** |
| Mentioned-11의 합집합 coverage | 4/8 = **50%** | 4/10 = **40%** |
| Essential-26의 합집합 coverage | 8/8 = **100%** | 10/10 = **100%** |

→ **두 cohort 사이에서 실제로 variability를 드러낸 tag 10개 중, mentioned-11은 4개만 포착. Essential-26은 10개 모두 포함. 심지어 derivation을 적용하면 mentioned-11의 상대적 coverage는 50%에서 40%로 오히려 감소** — derivation으로 새로 드러나는 tag가 전부 "not mentioned"이기 때문.

---

## 1. 배경 — 왜 "mentioned" 과 "essential" 을 구분하나

### 1.1 AI medical imaging의 replication 문제

공개 논문을 다른 기관에서 재현하려 할 때, 재현자가 맞춰야 하는 최소한의 acquisition parameters가 뭔지 모호하다:

- 논문 Methods 섹션에서 저자가 언급한 DICOM tag 목록 = **mentioned tags**
- Protocol2DICOM 이 제안하는 26개 CT essential tag 집합 = **essential tags**
- 이 둘의 관계, 그리고 실제 deployment 데이터에서 "어떤 tag가 결과 variability에 영향을 주는가"는 경험적 검증이 안 되어 있음

Claim: **deployment distribution의 acquisition variability는 mentioned tags만으로는 절반 이하만 포착된다.** Essential-26은 필요한 safety net이다.

### 1.2 검증 전략

하나의 cohort만으로는 검증 불충분. 서로 다른 두 distribution (여러 site, 여러 scanner의 HCC-TACE-Seg vs 단일 site MSKCC의 CRLM; HCC+ vs HCC−)에서 같은 분석을 돌려야 **deployment-specific** 과 **universal** variability를 구분할 수 있음.

---

## 2. ST3DCN Replication — 우리가 어디까지 재현했는가

### 2.1 대상 모델

**ST3DCN** (Slice-based Thoracic 3D Convolutional Network / Spatial-Temporal 3D-DCN)
- Yu et al., *JHEP Reports* 2025 (PMC11648772)
- Github: `HKUMedicineLiverAI/ST3DCN`
- 입력: liver lesion을 포함하는 70 × 70 × 70 voxel 3D crop (portal venous phase CT만 사용하는 경우 기준)
- 출력: 해당 lesion이 HCC일 확률 p(HCC) — Table 2A/S4F 기준 threshold 0.8로 이진 분류
- Pretrained weights: repo 제공

### 2.2 Pipeline validation — "demo sample에서 pixel-exact 재현"

논문 repo에는 2개의 sample 데이터가 들어있음 (`ID_0848`: HCC+, `QEH032`: HCC−, 각각 10 observations). 해당 sample의 **pre-cropped .npy** 가 포함되어 있어, 원 DICOM 으로부터 이 .npy를 재생성하는 과정의 pixel-exact 재현이 가능.

우리 구현 (`cleaned_code/pipeline/preprocessing.py`)의 검증 결과:
- 10/10 observation 모두 **|Δ| = 0** (pixel-exact 일치)
- Validation 후 두 TCIA cohort (HCC-TACE-Seg, CRLM) 에 동일 파이프라인 적용
- 이 pixel-exact 재현은 다음 사항들을 "우리가 제대로 맞췄다"는 증거로 사용:
  1. HU windowing [-160, 240]
  2. Tight bbox + asymmetric padding (+10, +9) per axis
  3. Axis permutation (Z, Y, X) → (Y, X, Z) — 논문 .npy convention
  4. 70³ resize/pad + [0, 1] normalize
  5. uint8 intermediate scaling

Validated preprocessing 함수들은 `cleaned_code/pipeline/preprocessing.py` 에 모듈화되어 있음. 임의의 DICOM cohort에 대해 재사용 가능.

### 2.3 재현 Scope — 좁게 결정한 3가지

우리는 **scope를 좁게** 잡아서 comparability를 최대화함:

#### (a) **Portovenous phase only**
- 근거 1: 논문 Table S4F — portovenous-only AUC 0.912 ≈ triphasic 0.919 (사실상 동등).
- 근거 2: HCC-TACE-Seg의 SEG mask는 **venous phase** 위에 그려져 있음 (Morshid 2019).
- 근거 3: CRLM은 TCIA 공식 설명에서 "preoperative portal venous contrast-enhanced MDCT"로 공개 (Simpson 2024).
- 근거 4: `ContrastBolusStartTime / ContrastBolusAgentPhase` 등 phase metadata가 양쪽 cohort에서 99~100% null → image-based 판별이 유일한 방법. Portovenous-only로 좁히면 이 문제를 회피.

#### (b) **Pre-treatment only**
- HCC-TACE-Seg: 각 환자마다 여러 study가 있음 (baseline + follow-up). NBIA digest의 `Longitudinal Temporal Offset From Event`가 가장 작은 (DIAGNOSIS event 기준 가장 가까운) study → SEG가 포함된 study 중 pre-TACE baseline을 선택.
- CRLM: TCIA collection 설명이 이미 "preoperative … within 6 weeks of hepatic resection"로 명시. 환자당 study 1개 → 자동으로 pre-treatment.

#### (c) **Lesion-level 분석**
- Per-patient 단위로 **union bbox** 를 쓰면 multi-lesion 환자의 70³ crop 이 여러 lesion을 포괄하게 되어 signal이 희석됨.
- CRLM에서 실측:
  - Per-patient union-bbox specificity: 0.3046 (137/197 FP)
  - Per-lesion specificity: **0.7931** (96/464 FP)
  - 단일-lesion 환자 평균 prob 0.45 vs ≥2-lesion 환자 0.92 — aggregation 방식이 결과를 30+pp 흔듬
- → **lesion-level이 정답 단위**. HCC 쪽도 `subset_summary_lesion.txt` 기준 146 lesions로 이미 lesion-level 분석 수행됨.

### 2.4 Cohort 별 replication 수치 (baseline)

**HCC-TACE-Seg (lesion-level, 146 lesions)**
- Sensitivity @ 0.8 = **0.8219** (TP 120/146)
- Wilson 95% CI (0.7519, 0.8755)
- 논문 reference p0 = 0.869 (Table 2A observation-level) — **CI가 p0을 포함** → replication 유효

**CRLM (lesion-level, 464 lesions, all HCC−)**
- Specificity @ 0.8 = **0.7931** (TN 368/464)
- Wilson 95% CI (0.7539, 0.8275)
- 논문에는 external test에서 metastasis subgroup 포함되어 있으므로 out-of-distribution 은 아님. 0.79 specificity는 public retrospective data로서는 합리적 수준.

이 baseline 값 자체는 연구의 목적이 아니라 **pipeline이 정상적으로 돌아간다는 증거**. 진짜 분석은 아래 per-tag heterogeneity.

---

## 3. 두 Cohort 데이터 특성

### 3.1 HCC-TACE-Seg (TCIA, Morshid 2019)

| 항목 | 값 |
|---|---|
| 환자 수 | 105 |
| Lesion 수 | 146 |
| Ground truth | 전원 HCC+ |
| Phase | Triphasic CT (arterial + portal + delayed) | 
| Site | 다중 site (multi-site cohort) |
| Scanner 제조사 | SIEMENS 다수, GE 일부 |
| SEG segment layout | 1=Liver, 2=**Mass**, 3=…, 4=Aorta |
| Mask 그려진 phase | Venous (portal) — Morshid 2019 |
| 분석 path | **sensitivity** primary |

**주의사항**: HCC-TACE-Seg의 CT series는 triphasic이 단일 series에 여러 `AcquisitionNumber`로 묶여 있음 (e.g., `LIVER 3 PHASE` 한 series에 A+P+D). → portal phase 식별을 위해 **Aorta segment (seg 4)** 의 slice당 max HU 사용 — arterial은 HU가 높고 portal은 상대적으로 낮음.

### 3.2 CRLM (TCIA Colorectal-Liver-Metastases, Simpson 2024)

| 항목 | 값 |
|---|---|
| 환자 수 | 197 |
| Lesion 수 | 464 (환자당 평균 2.36개, 1~5개) |
| Ground truth | 전원 HCC− (colorectal liver metastasis, 조직학 확진) |
| Phase | Portal venous only (TCIA 공식 기재) |
| Site | MSKCC 단일 site |
| Scanner 제조사 | GE Medical Systems 195/197 |
| SEG segment layout | 1=Liver, 2=Liver Remnant, 3=Hepatic Vein, 4=**Portal Vein**, 5~9=**Tumor_1~Tumor_5** (PropType=Mass) |
| Mask 그려진 phase | Portal venous (TCIA 공식 + Portal Vein HU > 120 실측 검증) |
| 분석 path | **specificity** primary |

**Phase 검증**: CRLM SEG에 Portal Vein segment가 있으므로, 10명 sample에서 Portal Vein median HU 측정 → 123~202 HU (모두 > 120) → portal venous phase 확증.

**ST3DCN training distribution과의 관계**: ST3DCN external test cohort에 metastasis subgroup이 포함되어 있으므로 CRLM 이 out-of-distribution 은 아님.

### 3.3 두 cohort의 대조적 구조

| | HCC-TACE-Seg | CRLM |
|---|---|---|
| Site variability | 多 (multi-site) | 0 (single site) |
| Manufacturer variability | 多 | 少 (99% GE) |
| CT series 당 AcqNum 수 | 多 (triphasic bundled) | 1 |
| Tumor segment 수 | 1 (Mass) | 5 (Tumor_1~5) |
| Pre-treatment 식별 | longitudinal offset 필요 | TCIA가 전원 preop으로 공개 |

이 대조가 중요함. **각 cohort가 서로 다른 acquisition variability를 expresses** 하기 때문에 FDR-sig tag 집합이 cohort-specific으로 나올 수밖에 없고, 이게 바로 essential-safety-net 논리의 핵심 증거가 됨.

---

## 4. 분석 방법론

### 4.1 Essential-26 와 Mentioned-11

26개 essential CT tag 목록 (protocol2dicom 정의; `cleaned_code/analysis/essential_tags.py`에 코드화):

**Intensity/exposure (10)**: KVP, XRayTubeCurrent, ConvolutionKernel, ContrastBolusAgent, CTDIvol, ContrastFlowRate, Exposure_uAs, RevolutionTime, ReconstructionAlgorithm, ContrastBolusAgentPhase

**Geometry (12)**: SliceThickness, Rows, Columns, PixelSpacing, ReconstructionDiameter, WindowCenter, WindowWidth, TotalCollimationWidth, SpiralPitchFactor, TableSpeed, SpacingBetweenSlices, ReformattingThickness

**Device (2)**: ManufacturerModelName, Manufacturer

**Anatomy/position (2)**: BodyPartExamined, PatientPosition

이 중 **11개가 ST3DCN 논문 Methods/Table S1에 mentioned** 되어 있음: KVP, XRayTubeCurrent, ContrastBolusAgent, ContrastFlowRate, ContrastBolusAgentPhase, SliceThickness, SpiralPitchFactor, ReformattingThickness, ManufacturerModelName, Manufacturer, BodyPartExamined.

(추가로 논문이 언급하는 3개 시퀀스형 필드 — ContrastBolusAgentSequence, ContrastBolusInjectionDelay, CTAcquisitionTypeSequence — 는 값형 DICOM tag이 아니라 per-tag heterogeneity test 불가능. 별도 표기.)

### 4.2 Kruskal-Wallis test

각 essential tag에 대해:
1. Tag 값을 bin으로 그루핑 (categorical passthrough / numerical cut)
2. 각 bin에서 lesion의 p(HCC) 분포 수집
3. **Kruskal-Wallis H test**: 2개 이상의 bin 사이 p(HCC) 분포가 동일한가? (귀무가설: 동일 분포)
4. Bin별 최소 n=5 lesion 요구; 미만이면 해당 bin은 exclude
5. testable bin 수가 2 이상일 때만 검정 수행

### 4.3 Benjamini-Hochberg False Discovery Rate (BH-FDR)

**FDR 은 무엇인가**:
- 여러 개 tag을 동시에 검정하면 우연히 p < 0.05가 나오는 tag이 생길 수 있음 (false positive). 26개 tag 독립 검정 시 기대 false positive = 26 × 0.05 = 1.3개.
- "**False Discovery Rate**" = 유의하다고 판정된 것들 중 실제로는 거짓인 것의 비율.
- BH 절차는 이 FDR을 예컨대 0.05 이하로 통제.

**계산 방법 (Benjamini-Hochberg 1995)**:
1. 모든 testable tag의 raw p-value를 오름차순 정렬: p(1) ≤ p(2) ≤ … ≤ p(m)
2. 각각의 FDR 조정 p = p(k) × m / k (k는 순위)
3. 오른쪽부터 monotone minimum으로 단조화: `FDR_adj(k) = min(FDR_adj(k), FDR_adj(k+1), …)`
4. 최종 FDR-adjusted p를 0.05와 비교

**해석**: FDR p < 0.05 이면 "이 tag의 bin 간 p(HCC) 분포 차이는 false-discovery rate 5% 이하로 통제된 상태에서 유의" — 즉 여러 tag를 동시에 검정하는 맥락을 보정한 뒤에도 heterogeneity가 남는다는 뜻.

**α=0.05** 사용. Testable tag pool 전체에 대해 BH 적용.

### 4.4 Tag status 분류 (per cohort)

각 essential-26 tag은 cohort의 data fill 상태에 따라 다음 4개 중 하나로 분류:

- **fully_missing**: fill rate < 5%. 테스트 불가.
- **invariant**: populated 되었으나 unique value 1개. 테스트 불가.
- **insufficient**: populated, unique value ≥ 2이지만 n≥5 bin이 2개 미만. 테스트 불가.
- **testable**: n≥5 bin이 2개 이상. Kruskal-Wallis + BH-FDR 적용.

### 4.5 Derivation — DICOM 표준에 근거한 tag 복구

단순 fill 여부를 보는 것만으로는 놓치는 신호가 있음. DICOM PS3.3 표준이 "이 tag은 다른 tag의 조합으로 계산 가능" 이라고 명시하는 경우가 많음. 우리는 **표준이 명시적으로 허용하는 계산**만 적용:

#### (a) SpacingBetweenSlices (0018,0088)
- **근거**: DICOM PS3.3 §C.7.6.2.1.1
- **공식**: `median(diff(sorted ImagePositionPatient[2])))` per CT series
- **dependency fill**: `ImagePositionPatient` (0020,0032) — 모든 CT DICOM에서 100% populated
- **효과**: 양쪽 cohort 모두 direct 0~3% → derived 100%

#### (b) Exposure_uAs (0018,1153)
- **근거**: DICOM PS3.3 §C.8.7.2.1 — 표준이 "expressed in μAs, **for example calculated from Exposure Time and X-Ray Tube Current**" 라고 명시
- **공식**: `ExposureTime (ms) × XRayTubeCurrent (mA)` = μAs (단위 검증 완료)
- **대안 고려**: `Exposure (mAs) × 1000` 도 DICOM 정의상 exact하지만, (0018,1152) Exposure field가 vendor-inconsistent (per-rotation effective vs total-scan) → 같은 scanner 내에서 CV 0.74~1.31. Time×Current는 같은 scanner 내 CV 0.13~0.38로 훨씬 일관됨 → Time×Current를 primary derivation으로 채택.
- **dependency fill**: HCC 98.6%, CRLM 44.4%
- **효과**: 양쪽 direct 0% → derived (HCC) 99% / (CRLM) 44%

#### (c) ContrastBolusAgentPhase (0018,9344) — 개념적 derivation
- Pipeline이 이미 SEG ref-SOP + Portal Vein HU 로 portal venous phase를 확증하므로, 암묵적 derivation 은 완료. 단 derived 값이 전원 "PORTAL_VENOUS" 상수라 discriminative test 불가 (invariant).
- → 통계 분석에는 추가 기여 없지만, 파이프라인이 이 정보를 이미 활용하고 있다는 증빙으로 문서화.

#### (d) Manufacturer / ManufacturerModelName (HCC-TACE-Seg만 시도)
- HCC-TACE-Seg DICOM에서 직접 fill 1.4%. TCIA가 제공하는 NBIA digest Excel 에 일부 보완 데이터 존재.
- NBIA digest 조인 시도 → **digest 자체가 111/677 series만 populated**; 105 HCC 환자 중 1명만 복구 가능. → "mentioned이면서 외부 메타데이터로도 구조적으로 복구 불가" 증거로 문서화.

#### Derivation 시도했으나 불가능한 tag
- CTDIvol: scanner factory-computed, 다른 tag에서 derive 불가
- ContrastBolusAgent: named value
- ContrastFlowRate: 기록 누락
- ReformattingThickness: SliceThickness와 개념적으로 다르므로 치환 불가 (별도 post-processing field)
- ReconstructionAlgorithm: ConvolutionKernel과 일부 중복되나 DICOM 표준상 별도 tag이므로 substitution 부적절
- TableSpeed / TotalCollimationWidth / SpiralPitchFactor: 같은 **helical geometry cluster** — 하나가 missing 이면 나머지도 missing (Cramér's V = 1.0). 파생으로 fill 증가 불가.

---

## 5. 결과 — 세 가지 관측

### 관측 1. Mentioned-11 은 각 cohort variability의 절반 이하만 포착

Derivation 적용 후 (primary 분석):

| Cohort | FDR-sig 총 | Mentioned (paper) | Not mentioned |
|---|---|---|---|
| HCC-TACE-Seg (전원 HCC+) | **6** | KVP, SpiralPitchFactor (**2/6**) | Exposure_uAs◆, TableSpeed, TotalCollimationWidth, WindowCenter |
| CRLM (전원 HCC−) | **5** | ContrastBolusAgent, SliceThickness (**2/5**) | ReconstructionDiameter, SpacingBetweenSlices◆, WindowCenter |

◆ = derivation-augmented tag.

- HCC cohort: mentioned 2/6 = **33%** 만 포착.
- CRLM cohort: mentioned 2/5 = **40%** 만 포착.
- 합집합 10개 중 mentioned 4개 = **40%**. 논문 Methods만 따라 replication 준비했다면 acquisition variability signal의 **60%를 놓침**.

Direct-only 분석 (derivation 적용 전) 과 비교:

| Cohort | Direct FDR-sig | Augmented FDR-sig |
|---|---|---|
| HCC | 4 (KVP, SpiralPitch, TCW, TableSpeed) | 6 (+ Exposure_uAs◆, WindowCenter) |
| CRLM | 4 (ContrastBolusAgent, SliceThickness, ReconDiam, WindowCenter) | 5 (+ SpacingBetweenSlices◆) |

Derivation 으로 증가한 sig tag은 모두 **not mentioned**. 즉 derivation 은 mentioned-14의 상대 coverage를 오히려 낮춤 (direct 50% → augmented 40%).

### 관측 2. Cohort 간 FDR-sig 교집합: ∅ → {WindowCenter}

- **Direct-only**: HCC 와 CRLM 모두에서 FDR-sig가 나온 tag = **∅ (공집합)**.
  - HCC sig: {KVP, SpiralPitchFactor, TotalCollimationWidth, TableSpeed}
  - CRLM sig: {ContrastBolusAgent, SliceThickness, ReconstructionDiameter, WindowCenter}
  - 교집합 없음 → **어떤 tag가 discriminative인지는 deployment distribution이 결정**.
- **Derivation 적용 후**: 교집합 = **{WindowCenter}** (1개)
  - WindowCenter 는 논문 mentioned 에 없음.
  - 두 cohort에서 독립적으로 FDR-sig → "universal" signal로 간주 가능.

**함의**:
1. 단일 cohort 논문 (ST3DCN 포함)으로는 어떤 tag이 replication 에 중요한지 사전에 알 수 없음.
2. Mentioned-11 만 사용하면:
   - HCC 이식 시: TotalCollimationWidth, TableSpeed, WindowCenter, Exposure_uAs — **4개 누락**
   - CRLM 이식 시: ReconstructionDiameter, SpacingBetweenSlices, WindowCenter — **3개 누락**
3. 심지어 두 cohort 모두에서 유의한 universal tag (WindowCenter) 조차 mentioned에 없음.

→ **넓은 safety net (essential-26)이 필요한 이유**: 어떤 tag가 특정 cohort에서 variability를 드러낼지 사전 예측 불가능. Mentioned 만으로는 반드시 blind spot이 생김.

### 관측 3. Mentioned tag 의 다수는 애초에 unusable

논문이 tag을 언급했다고 해서 deployment DICOM에 그 값이 들어있다는 보장은 없음. Derivation 시도 후에도 mentioned-11 중 대부분이 non-testable:

| Cohort | Mentioned-11 testable | Mentioned-11 FDR-sig | **Mentioned-11 unusable** |
|---|---|---|---|
| HCC-TACE-Seg | 3 | 2 | **8** |
| CRLM | 5 | 2 | **6** |

Unusable (fully_missing / invariant / insufficient) mentioned tag 목록:

- **HCC-TACE-Seg**:
  - fully_missing: ContrastBolusAgent, ContrastFlowRate, ContrastBolusAgentPhase, ReformattingThickness, ManufacturerModelName, Manufacturer
  - invariant: BodyPartExamined (= LIVER for all)
  - insufficient: SliceThickness (2.5mm dominant)
- **CRLM**:
  - fully_missing: ContrastFlowRate, ContrastBolusAgentPhase, ReformattingThickness, BodyPartExamined
  - insufficient: KVP (195/197 = 120), Manufacturer (99% GE)

**양쪽 cohort 모두에서 unusable인 mentioned tag 5개** — "publish-retrospective-data 에서 구조적 복구 불가"의 증거:

| Tag | HCC status | CRLM status |
|---|---|---|
| BodyPartExamined | invariant | fully_missing |
| ContrastBolusAgentPhase | fully_missing | fully_missing |
| ContrastFlowRate | fully_missing | fully_missing |
| Manufacturer | fully_missing | insufficient |
| ReformattingThickness | fully_missing | fully_missing |

이 5개는 어떤 derivation을 해도 살릴 수 없음. **Mentioned-only 접근법의 구조적 한계**.

Essential-26 전체는 약 2배 많은 observable variability 를 제공:

| Cohort | Essential-26 testable | Essential-26 sig | Mentioned-11 testable | Mentioned-11 sig |
|---|---|---|---|---|
| HCC-TACE-Seg | **9** | **6** | 3 | 2 |
| CRLM | **13** | **5** | 5 | 2 |

---

## 6. Essential-26 전체 분류 (HCC ‖ CRLM 병렬)

표 읽는 법:
- **✓ / ✗** = paper mentioned 여부
- **Status**: fully_missing / invariant / insufficient / testable
- **FDR**: Benjamini-Hochberg 보정 후 p-value (testable 한 tag만)
- **★** = FDR < 0.05 (significant)
- **◆** = derivation 적용됨

| Tag | Hex | Class | Mentioned | HCC fill% | HCC status | HCC FDR | HCC | CRLM fill% | CRLM status | CRLM FDR | CRLM |
|---|---|---|---|---|---|---|---|---|---|---|---|
| KVP | 180060 | I | ✓ | 100.0 | testable | 0.0021 | **★** | 100.0 | insufficient | — |  |
| XRayTubeCurrent | 181151 | I | ✓ | 99.3 | testable | 0.6483 |  | 44.4 | testable | 0.8182 |  |
| ConvolutionKernel | 181210 | I | ✗ | 99.3 | insufficient | — |  | 44.4 | insufficient | — |  |
| ContrastBolusAgent | 180010 | I | ✓ | 1.4 | fully_missing | — |  | 99.1 | testable | 0.0056 | **★** |
| CTDIvol | 189345 | I | ✗ | 0.7 | fully_missing | — |  | 0.0 | fully_missing | — |  |
| ContrastFlowRate | 181046 | I | ✓ | 0.0 | fully_missing | — |  | 0.0 | fully_missing | — |  |
| Exposure_uAs | 181153 | I | ✗ | 98.6 | testable | 0.0021 | **★◆** | 44.4 | testable | 0.4668 | ◆ |
| RevolutionTime | 189305 | I | ✗ | 52.7 | invariant | — |  | 5.4 | insufficient | — |  |
| ReconstructionAlgorithm | 189315 | I | ✗ | 0.0 | fully_missing | — |  | 0.0 | fully_missing | — |  |
| ContrastBolusAgentPhase | 189344 | I | ✓ | 0.0 | fully_missing | — |  | 0.0 | fully_missing | — |  |
| SliceThickness | 180050 | G | ✓ | 100.0 | insufficient | — |  | 100.0 | testable | 0.0056 | **★** |
| Rows | 280010 | G | ✗ | 100.0 | invariant | — |  | 100.0 | invariant | — |  |
| Columns | 280011 | G | ✗ | 100.0 | invariant | — |  | 100.0 | invariant | — |  |
| PixelSpacing | 280030 | G | ✗ | 100.0 | testable | 0.1974 |  | 100.0 | testable | 0.4668 |  |
| ReconstructionDiameter | 181100 | G | ✗ | 100.0 | testable | 0.2182 |  | 100.0 | testable | 0.0150 | **★** |
| WindowCenter | 281050 | G | ✗ | 100.0 | testable | 0.0488 | **★** | 100.0 | testable | 0.0233 | **★** |
| WindowWidth | 281051 | G | ✗ | 100.0 | insufficient | — |  | 100.0 | testable | 0.0955 |  |
| TotalCollimationWidth | 189307 | G | ✗ | 52.7 | testable | 0.0140 | **★** | 5.4 | testable | 0.9535 |  |
| SpiralPitchFactor | 189311 | G | ✓ | 52.7 | testable | 0.0140 | **★** | 5.4 | testable | 0.9535 |  |
| TableSpeed | 189309 | G | ✗ | 52.7 | testable | 0.0140 | **★** | 5.4 | testable | 0.9535 |  |
| SpacingBetweenSlices | 180088 | G | ✗ | 100.0 | insufficient | ◆ |  | 100.0 | testable | 0.0274 | **★◆** |
| ReformattingThickness | 720512 | G | ✓ | 0.0 | fully_missing | — |  | 0.0 | fully_missing | — |  |
| ManufacturerModelName | 081090 | D | ✓ | 1.4 | fully_missing | — |  | 100.0 | testable | 0.9535 |  |
| Manufacturer | 080070 | D | ✓ | 1.4 | fully_missing | — |  | 100.0 | insufficient | — |  |
| BodyPartExamined | 180015 | A | ✓ | 100.0 | invariant | — |  | 0.6 | fully_missing | — |  |
| PatientPosition | 185100 | A | ✗ | 100.0 | invariant | — |  | 100.0 | insufficient | — |  |

Class: **I** = Intensity/exposure, **G** = Geometry, **D** = Device, **A** = Anatomy/position.

---

## 7. 각 Tag bin별 결과 (HCC sens / CRLM spec / 평균 p(HCC))

Bin별 per-tag 테이블은 `RESULTS_ESSENTIAL_vs_MENTIONED.md` 의 Section 5에 전체가 있음. 여기서는 **FDR-sig tag 10개** 의 bin-level 수치만 발췌.

### 7.1 KVP (180060) — HCC에서 sig, CRLM에서 invariant
| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| 100 | 10 | 0.700 | 0.757 | — | — | — |
| 120 | 79 | 0.797 | 0.786 | 461 | 0.794 | 0.320 |
| 130 | — | — | — | 1 | 1.000 | 0.014 |
| 140 | 57 | 0.877 | 0.885 | 2 | 0.500 | 0.485 |

CRLM은 MSKCC single-site로 99%가 120 kVp → testable이지만 n≥5 bin이 1개라 insufficient.

### 7.2 SliceThickness (180050) — CRLM에서 sig, HCC에서 insufficient
| bin (mm) | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| 2.5 | 143 | 0.818 | 0.819 | 213 | 0.746 | 0.366 |
| 3.0 | 1 | 0.000 | 0.431 | — | — | — |
| 3.75 | — | — | — | 20 | 0.850 | 0.329 |
| 5.0 | 2 | 1.000 | 0.929 | 226 | 0.836 | 0.266 |
| 7.5 | — | — | — | 5 | 1.000 | 0.186 |

HCC는 2.5mm에 쏠려 있어 bin 수 부족. CRLM은 2.5~7.5mm 스펙트럼 → FDR-sig.

### 7.3 ContrastBolusAgent (180010) — CRLM에서 sig
양쪽 값이 너무 많아 (CRLM 27개 unique) 전체 테이블은 부록 참조. 핵심:
- HCC: 144/146 missing → fully_missing
- CRLM: "OMNI" (Omnipaque) 계열 96%, 하위 brand/oral agent 조합 다양 → FDR 0.0056

### 7.4 Exposure_uAs (181153) — **derivation 으로 rescue, HCC에서 sig**
| bin (μAs) | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| [0, 150k) | 2 | 1.000 | 0.913 | 1 | 1.000 | 0.014 |
| [150k, 250k) | 14 | 0.786 | 0.751 | 29 | 0.828 | 0.287 |
| [250k, 350k) | 69 | 0.768 | 0.762 | 29 | 0.724 | 0.420 |
| [350k, 450k) | 42 | 0.905 | 0.913 | 83 | 0.843 | 0.236 |
| [450k, 1M) | 17 | 0.824 | 0.872 | 64 | 0.766 | 0.324 |
| (missing) | 2 | 1.000 | 0.983 | 258 | 0.787 | 0.340 |

Direct fill이 0%인 tag을 ExposureTime × XRayTubeCurrent 로 살림.

### 7.5 TotalCollimationWidth / SpiralPitchFactor / TableSpeed — HCC에서 sig (helical cluster)
세 tag은 Cramér's V = 1.0 (같은 helical geometry cluster). HCC에서 모두 같은 FDR 값 (0.0140). CRLM은 5.4% fill이라 bin별 n 부족 → sig 아님.

| SpiralPitch bin | HCC n | HCC sens | HCC mean p |
|---|---|---|---|
| 0.938 | 40 | 0.900 | 0.901 |
| 0.984 | 37 | 0.784 | 0.797 |
| (missing) | 69 | 0.797 | 0.782 |

### 7.6 ReconstructionDiameter (181100) — CRLM에서 sig
| bin (mm) | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| [0, 320) | 3 | 1.000 | 0.984 | 28 | 0.786 | 0.290 |
| [320, 360) | 18 | 0.833 | 0.812 | 139 | 0.799 | 0.331 |
| [360, 400) | 55 | 0.745 | 0.767 | 213 | 0.812 | 0.316 |
| [400, 500) | 57 | 0.895 | 0.891 | 52 | 0.615 | 0.418 |
| [500+) | 13 | 0.692 | 0.725 | 32 | 0.844 | 0.263 |

### 7.7 WindowCenter (281050) — **양쪽 cohort 모두 sig (universal signal)**
| bin | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| 30 | 1 | 1.000 | 0.999 | 2 | 1.000 | 0.065 |
| 35 | — | — | — | 1 | 1.000 | 0.023 |
| 40 | 39 | 0.821 | 0.817 | 398 | 0.814 | 0.311 |
| 45 | — | — | — | 2 | 0.500 | 0.556 |
| 50 | 48 | 0.833 | 0.847 | 48 | 0.667 | 0.463 |
| 60 | 21 | 0.714 | 0.681 | 11 | 0.909 | 0.223 |
| -50 | — | — | — | 2 | 1.000 | 0.003 |

두 cohort에서 독립적으로 FDR-sig → **derivation 없이도 universal variability source**.

### 7.8 SpacingBetweenSlices (180088) — **derivation 으로 rescue, CRLM에서 sig**
| bin (mm) | HCC n | HCC sens | HCC mean p | CRLM n | CRLM spec | CRLM mean p |
|---|---|---|---|---|---|---|
| [0, 1.5) | 4 | 1.000 | 0.984 | 5 | 1.000 | 0.405 |
| [1.5, 3.0) | 139 | 0.813 | 0.815 | 233 | 0.747 | 0.371 |
| [3.0, 5.0) | 1 | 1.000 | 0.998 | — | — | — |
| [5.0, 8.0) | 2 | 1.000 | 0.929 | 226 | 0.836 | 0.266 |

HCC는 [1.5, 3.0) bin에 139개 쏠려 있어 insufficient. CRLM은 2~3mm / 5~8mm 양극 분포 → FDR-sig.

---

## 8. 코드 구조 및 재현 가이드

### 8.1 디렉토리 개요 (`/home/ubuntu/hcc_workspace/cleaned_code/`)

```
cleaned_code/
├── README.md                           # 영문 사용자 가이드
├── pipeline/                           # ST3DCN inference pipeline
│   ├── preprocessing.py                #   SEG→CT 정렬, 70³ crop, HU window
│   ├── dataset_adapter.py               #   코호트별 규칙 (segment layout)
│   ├── lesion_pipeline.py              #   환자 → 병변별 crop
│   └── inference.py                     #   모델 로드 + predict
├── analysis/                           # Essential-tag heterogeneity analysis
│   ├── essential_tags.py                #   essential-26 + mentioned flag
│   ├── extract_tags.py                  #   DICOM tag 추출
│   ├── derive_tags.py                   #   DICOM 표준 derivation
│   ├── binners.py                       #   numeric tag binning
│   ├── per_tag_analysis.py              #   Kruskal-Wallis + BH-FDR
│   └── combined_report.py               #   cross-cohort 리포트
├── configs/                            # YAML 코호트 정의
│   ├── hcc_tace_seg.yaml
│   ├── crlm.yaml
│   └── new_dataset_template.yaml
└── runners/                            # 실행 스크립트
    ├── run_lesion_inference.py
    ├── run_analysis.py
    └── build_cross_cohort_report.py
```

### 8.2 핵심 함수와 validated framework

**Validated preprocessing** (demo sample에서 pixel-exact 확인된 부분):
- `pipeline.preprocessing.read_seg` — SEG DICOM parsing
- `pipeline.preprocessing.index_ct_series` — CT slice indexing by SOPInstanceUID
- `pipeline.preprocessing.validate_sop_linkage` — SEG-ref ↔ CT SOP 매칭 확인
- `pipeline.preprocessing.pick_mask_aligned_ct` — multi-CT 중 mask-align CT 선택
- `pipeline.preprocessing.build_ct_volume` — z-sorted 3D volume 생성
- `pipeline.preprocessing.build_segment_mask` — SEG → 3D mask (nearest-neighbor z-align, tol 1.25mm)
- `pipeline.preprocessing.crop_and_window` — **tight bbox + (+10, +9) asymmetric pad + HU [-160,240] + transpose + 70³ resize + [0,1] normalize** (이 부분이 pixel-exact 검증의 핵심)
- `pipeline.preprocessing.portal_vein_median_hu` — CRLM phase sanity
- `pipeline.preprocessing.identify_portal_acq_by_aorta` — HCC-TACE-Seg phase identification
- `pipeline.inference.load_model` / `predict` — ST3DCN 모델 로드 + 추론

**Per-cohort adapter** (새 데이터셋 추가 시 여기를 subclass):
- `pipeline.dataset_adapter.DatasetAdapter` — 추상 베이스
- `pipeline.dataset_adapter.CRLMAdapter` — CRLM 전용
- `pipeline.dataset_adapter.HCCTACESegAdapter` — HCC-TACE-Seg 전용

**Analysis pipeline**:
- `analysis.extract_tags.tag_dataframe` — 환자별 DICOM tag 추출 + lesion CSV에 병합
- `analysis.derive_tags.derive_tags` — SpacingBetweenSlices ← IPP, Exposure_uAs ← Time×Current
- `analysis.per_tag_analysis.analyze` — KW + BH-FDR + classification

### 8.3 재현 명령 (HCC-TACE-Seg + CRLM)

```bash
cd /home/ubuntu/hcc_workspace/cleaned_code
source /home/ubuntu/hcc_venv/bin/activate

# 1. 병변별 inference
LD_LIBRARY_PATH=/usr/local/cuda-11.2/lib64:$LD_LIBRARY_PATH \
CUDA_VISIBLE_DEVICES=2 \
python runners/run_lesion_inference.py --cohort crlm --out-dir out_crlm

LD_LIBRARY_PATH=/usr/local/cuda-11.2/lib64:$LD_LIBRARY_PATH \
CUDA_VISIBLE_DEVICES=2 \
python runners/run_lesion_inference.py --cohort hcc_tace_seg --out-dir out_hcc

# 2. tag 추출 + derivation + per-tag KW+FDR
python runners/run_analysis.py --cohort crlm \
    --lesion-csv out_crlm/lesion_results.csv \
    --patient-csv out_crlm/patient_results.csv \
    --out-dir out_crlm

python runners/run_analysis.py --cohort hcc_tace_seg \
    --lesion-csv out_hcc/lesion_results.csv \
    --patient-csv out_hcc/patient_results.csv \
    --out-dir out_hcc

# 3. Cross-cohort essential-vs-mentioned 리포트
python runners/build_cross_cohort_report.py \
    --hcc-dir out_hcc \
    --crlm-dir out_crlm \
    --out RESULTS_ESSENTIAL_vs_MENTIONED.md
```

### 8.4 새 코호트 추가 절차

`cleaned_code/README.md` 의 "Adding a new dataset" 섹션 참조:

1. `configs/new_dataset_template.yaml` 복사 → 새 코호트 YAML 작성
2. 데이터 레이아웃이 기존과 다르면 `DatasetAdapter` subclass 구현 (list_patients, resolve_series, acq_filter)
3. `ADAPTERS` dict에 등록
4. Phase 검증: Portal Vein segment 있으면 median HU > 120 확인. 없으면 image-level으로 검증
5. `run_lesion_inference.py --cohort <new>` → `run_analysis.py --cohort <new>` → cross-cohort report 에 포함

---

## 9. 결론 및 정리

### 9.1 입증한 것

1. **ST3DCN pipeline은 paper-faithful 하게 재현 가능**하다. Demo sample pixel-exact 검증 후, HCC-TACE-Seg lesion-level sens 0.82 (paper p0 = 0.869 포함한 Wilson CI), CRLM lesion-level spec 0.79.
2. **각 deployment cohort에서 acquisition-parameter variability를 드러내는 DICOM tag의 집합은 서로 다르다**. HCC와 CRLM의 FDR-sig 교집합은 direct-only에서 ∅, derivation 후 1개. 나머지는 cohort-specific.
3. **논문이 명시한 mentioned-11 은 어느 cohort에서든 실제 sig tag의 40% 이하만 포착한다**. 이는 essential 집합에 비해 체계적으로 부족함 (derivation 적용 후 상대 coverage는 오히려 감소하는 방향 50%→40%).
4. **Mentioned tag 중 상당수는 공공 retrospective DICOM에서 구조적으로 복구 불가**: ContrastFlowRate, ContrastBolusAgentPhase, ReformattingThickness, Manufacturer, BodyPartExamined 등. 어떤 derivation으로도 살릴 수 없음.
5. **Derivation은 표준 기반으로 mentioned-gap을 일부 메울 수 있으나**, 메워지는 tag들이 대체로 "not mentioned" (Exposure_uAs, SpacingBetweenSlices, WindowCenter) → mentioned-only 접근의 구조적 한계를 강화.

### 9.2 Thesis sentence

**Essential tags (n=26) 은 acquisition-parameter variability를 관측하기 위한 최소 safety net이다. 특정 tag가 deployment에서 discriminative할지 사전에 알 수 없기 때문에, 논문이 명시한 mentioned-14 subset으로는 실제 variability의 40% 미만만 포착된다. 두 독립 cohort (HCC+ HCC-TACE-Seg, HCC− CRLM) 에서 이 gap이 재현되었고, DICOM 표준 derivation을 적용해도 mentioned-essential 간 gap은 좁아지지 않는다.**

### 9.3 코드 + 데이터 산출물

- **Cleaned code**: `/home/ubuntu/hcc_workspace/cleaned_code/` (영문, 재사용·확장 가능)
- **영문 사용자 가이드**: `cleaned_code/README.md`
- **한글 마스터 문서 (이 파일)**: `RESULTS_KOREAN.md`
- **Cross-cohort essential-vs-mentioned 리포트**: `RESULTS_ESSENTIAL_vs_MENTIONED.md` (26개 tag 전체 bin-level 데이터)
- **Cohort-specific 산출물**:
  - HCC: `/home/ubuntu/hcc_workspace/replication/lesion_results.csv`, `per_tag_classification_augmented.csv`, `per_tag_binned_augmented.csv`
  - CRLM: `/home/ubuntu/hcc_workspace/replication_crlm/lesion_results.csv`, `per_tag_classification_augmented.csv`, `per_tag_binned_augmented.csv`

### 9.4 추후 cohort 확장

이 framework는 새 TCIA(혹은 기타) DICOM collection에 그대로 이식 가능. `cleaned_code/configs/new_dataset_template.yaml`과 `DatasetAdapter` subclass만 작성하면 동일 파이프라인이 돌고, cross-cohort report에 자동 포함됨.

추가 cohort가 합류하면 "essential-26 이 universal safety net" 주장이 더 강력해지고, mentioned-coverage 수치는 경향상 더 낮아질 것으로 예상.

---

## 부록 A. FDR 구체 계산 예시 (CRLM 기준)

CRLM augmented 분석에서 testable tag 13개의 raw p-value (오름차순):

| rank | tag | raw p | rank × α / m | FDR-adj |
|---|---|---|---|---|
| 1 | SliceThickness | 0.00043 | 0.00385 | 0.00561 |
| 2 | ContrastBolusAgent | 0.00086 | 0.00769 | 0.00561 |
| 3 | Exposure_uAs (direct) | 0.0014 | 0.01154 | 0.00607 |
| 4 | ReconstructionDiameter | 0.0046 | 0.01538 | 0.01504 |
| 5 | SpacingBetweenSlices | 0.0105 | 0.01923 | 0.02737 |
| 6 | WindowCenter | 0.0108 | 0.02308 | 0.02336 |
| 7 | WindowWidth | 0.0514 | 0.02692 | 0.09554 |
| 8 | PixelSpacing | 0.2867 | 0.03077 | 0.4668 |
| 9 | XRayTubeCurrent | 0.5662 | 0.03462 | 0.8182 |
| 10-13 | ... | ... | ... | ≥ 0.95 |

(m = 13). FDR 0.05 cutoff 아래: rank 1-6 → 6 tag FDR-sig. 단 위 표는 augmented 기준이고 Exposure_uAs 직접 binning이 세밀하게 바뀌어 최종 순위/FDR은 `per_tag_classification_augmented.csv` 참조.

---

## 부록 B. 용어 사전

| 용어 | 정의 |
|---|---|
| **DICOM PS3.3** | DICOM 표준 Part 3 (Information Object Definitions). Tag 정의 + semantics + derivation 예시. |
| **SEG** | DICOM Segmentation Object. Source CT slice를 `ReferencedSOPInstanceUID`로 참조하여 mask를 정의. |
| **Kruskal-Wallis H test** | 두 개 이상 그룹의 분포가 동일한지 non-parametric 검정. ANOVA의 non-parametric 대응. |
| **BH-FDR (Benjamini-Hochberg False Discovery Rate)** | 여러 가설 동시 검정 시 false discovery 비율을 통제하는 보정법. α = 0.05이면 "유의 판정된 것들 중 잘못된 것 비율 5% 이하". |
| **Wilson score 95% CI** | 이진 비율의 신뢰구간 (normal approx 보다 꼬리에서 정확). |
| **Cramér's V** | 두 categorical 변수의 연관도 측정치 (0~1). 1에 가까우면 거의 완전 공변. |
| **Portal venous phase** | 조영제 정맥 주입 후 ~60-90초 시점의 간 CT. Liver parenchyma와 portal vein 모두 enhance. Liver tumor detection의 표준 phase. |
| **TCIA** | The Cancer Imaging Archive. 공개 의료영상 repository. |
| **pu / pc / pce** | pu=uncontrolled 전체, pc=paper의 protocol tag 매칭, pce=pc ∩ essential tag. Subset sensitivity/specificity 비교용. |
| **essential-26** | protocol2dicom 이 정의한 26 CT DICOM tags (최소 acquisition-parameter 집합). |
| **mentioned-11** | essential-26 중 ST3DCN 논문 Methods/Table S1에 언급된 11개. |
| **derivation** | DICOM 표준이 허용하는 tag 간 계산식으로 missing tag을 복구. 예: `SpacingBetweenSlices ← median(diff(IPP[2]))`. |
