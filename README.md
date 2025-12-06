# IlluMeta

Illumina 메틸레이션 어레이(GEO 원시 IDAT) 다운로드 → 전처리 → minfi/sesame 분석 → HTML 대시보드 생성 파이프라인입니다.

## 초보자용 빠른 가이드 (Ubuntu 예시)
1) 필수 시스템 패키지
```bash
sudo apt-get update
sudo apt-get install -y python3 python3-pip r-base pandoc \
    libcurl4-openssl-dev libssl-dev libxml2-dev libicu-dev
```
2) Python 의존성 (IDAT 검색 스크립트용)
```bash
python3 -m pip install --upgrade pip requests
```
3) 저장소 위치에서 실행 (최초 1회 R 패키지 설치 포함)
```bash
# GEO 메타데이터/IDAT 내려받기
python3 illumeta.py download GSE12345 -o /path/to/project

# configure.tsv의 primary_group 컬럼을 수동으로 채움

# 분석 실행
python3 illumeta.py analysis -i /path/to/project --group_con Control --group_test Case
```
4) IDAT이 있는 GSE 검색 (옵션)
```bash
python3 illusearch.py --keywords "breast cancer" --email you@example.com -o geo_idat_methylation.tsv
```

## 필요 구성 요소(요약)
- Python 3 (illusearch.py는 반드시 `requests` 필요)
- R (`Rscript` 포함)
- pandoc (HTML 위젯 렌더링)
- 빌드용 시스템 라이브러리: `libcurl4-openssl-dev libssl-dev libxml2-dev` (Ubuntu)

## 자동 설치/캐시
- `illumeta.py` 최초 실행 시 `r_scripts/setup_env.R`가 자동 호출되어 필요한 R 패키지를 설치하고 sesame 데이터를 `cache/`에 캐시합니다. 한 번 성공하면 이후 실행은 설치 단계를 건너뜁니다.
- 강제로 다시 설치/검증하려면 실행 전에 `ILLUMETA_FORCE_SETUP=1`을 설정하세요.

## 자주 발생하는 설치 문제 해결
- `R is not installed or not in the PATH`: 위의 apt 명령으로 R 설치 후 터미널을 새로 열어 재실행.
- `library path not writable` 또는 R 패키지 설치 실패: 사용자 라이브러리 지정 후 재시도
  ```bash
  export R_LIBS_USER="$HOME/R/library"
  mkdir -p "$R_LIBS_USER"
  Rscript r_scripts/setup_env.R
  ```
- `No module named 'requests'`: `python3 -m pip install requests`
- `libxml2.so.2: cannot open shared object file` 등 xml2/XML 관련 오류:
  ```bash
  sudo apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev libicu-dev
  conda deactivate   # conda R이 간섭하면 비활성화
  export R_LIBS_USER="$HOME/R/library" && mkdir -p "$R_LIBS_USER"
  Rscript r_scripts/setup_env.R
  ```
- EPIC v2 패키지(매니페스트/어노테이션)가 필수입니다. Bioconductor/GitHub 이름(예: `IlluminaHumanMethylationEPICv2anno.20a1.hg38`, `IlluminaHumanMethylationEPICanno.ilm10b5.hg38`)을 자동 시도합니다. 실패하면 Bioc 업그레이드 또는 GitHub 소스 설치가 필요할 수 있습니다.

## GEO IDAT 검색 스크립트 (illusearch.py)
- 의존성: `python3 -m pip install requests`
- 사용 예시:
```bash
python3 illusearch.py --keywords "blood" --email you@example.com --retmax 500 -o geo_idat_methylation.tsv
```
`-o`를 생략하면 현재 디렉터리에 `geo_idat_methylation.tsv`로 저장됩니다.
- GEO 보조 파일(suppl)에 실제 IDAT/RAW가 있는지 기본적으로 확인합니다. 속도를 위해 생략하려면 `--no-check-suppl`을 사용하세요.
- TSV 컬럼: `gse_id`, `title`, `platform`(정규화된 GPL ID), `platform_type`(450k/850k/950k), `sample_count`, `pubmed_id`, `suppl_has_idat`(`yes`/`no`/`error`).

## 분석 파이프라인 결과물
- `*_results/` 폴더에 HTML/CSV 결과와 QC 메트릭스
- 상위 폴더에 대시보드 `Test_vs_Control_results_index.html`

## 주의 사항
- GEO 메타데이터에 원시 IDAT가 없으면 `download` 단계가 실패합니다. 이 경우 IDAT을 직접 `/path/to/project/idat/`에 배치하고 `configure.tsv`를 수동으로 준비한 뒤 `analysis`를 실행하세요.
- 캐시는 저장소 루트의 `cache/`에 저장됩니다. 다른 위치를 쓰려면 실행 시 현재 디렉터리를 변경해 실행하세요.
