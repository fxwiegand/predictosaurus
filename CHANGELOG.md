# Changelog

## [0.5.0](https://github.com/fxwiegand/predictosaurus/compare/v0.4.2...v0.5.0) (2025-09-02)


### Features

* Add --min-vaf parameter to build subcommand ([#141](https://github.com/fxwiegand/predictosaurus/issues/141)) ([5ca19c3](https://github.com/fxwiegand/predictosaurus/commit/5ca19c3a7c9cad5b03b4f68f3c062d25c887af62))


### Bug Fixes

* Add mutex to synchronize score writing in parallel loop ([#135](https://github.com/fxwiegand/predictosaurus/issues/135)) ([1f44ae2](https://github.com/fxwiegand/predictosaurus/commit/1f44ae226de3afbb71698b907d76751cc2c991be))
* Fix contains_variant to take target string argument ([#137](https://github.com/fxwiegand/predictosaurus/issues/137)) ([4a1a9e3](https://github.com/fxwiegand/predictosaurus/commit/4a1a9e3d02e48489a261e8c38cf0c8cbfa84f7eb))
* Gracefully handle ambiguous nucleotide 'N' for amino acid changes ([#133](https://github.com/fxwiegand/predictosaurus/issues/133)) ([f0e72a8](https://github.com/fxwiegand/predictosaurus/commit/f0e72a8cd5940e9bed4806014f2f0615199d67ff))
* Prevent cycles in forward traversal by checking path containment ([#136](https://github.com/fxwiegand/predictosaurus/issues/136)) ([31f0cf4](https://github.com/fxwiegand/predictosaurus/commit/31f0cf4c3e2b1f4e319cd3d1efc3599a6a2a976f))


### Performance Improvements

* Use transaction for writing graphs via DuckDB ([cdc2873](https://github.com/fxwiegand/predictosaurus/commit/cdc28730bd8d82eed60422b43c2139fce2ed8fa4))

## [0.4.2](https://github.com/fxwiegand/predictosaurus/compare/v0.4.1...v0.4.2) (2025-08-20)


### Bug Fixes

* Refactor get_targets to return only chromosomes that are actually ([#132](https://github.com/fxwiegand/predictosaurus/issues/132)) ([1e1c57b](https://github.com/fxwiegand/predictosaurus/commit/1e1c57b976a7e6c21f7708f226d72cefd36f2df6))


### Performance Improvements

* Add threads option and parallelize build and process command ([#130](https://github.com/fxwiegand/predictosaurus/issues/130)) ([8c83f96](https://github.com/fxwiegand/predictosaurus/commit/8c83f96c88e87616755ec46934a044c1db9dcd48))

## [0.4.1](https://github.com/fxwiegand/predictosaurus/compare/v0.4.0...v0.4.1) (2025-08-18)


### Bug Fixes

* Write score output to single TSV file instead of per transcript ([#124](https://github.com/fxwiegand/predictosaurus/issues/124)) ([df0773d](https://github.com/fxwiegand/predictosaurus/commit/df0773d3dfe8748754a4530deb2648f760216c08))

## [0.4.0](https://github.com/fxwiegand/predictosaurus/compare/v0.3.2...v0.4.0) (2025-08-06)


### Features

* Add amino acid distance metrics and Grantham matrix ([adfb9c1](https://github.com/fxwiegand/predictosaurus/commit/adfb9c11c3bc2dda45144c5673f51c3f09ca991c))
* Add distance matrix and calculation for sneath ([b09058c](https://github.com/fxwiegand/predictosaurus/commit/b09058cb43febfa33a272e627de428bfd6de578b))
* Add score storage and output for transcripts ([d0919f2](https://github.com/fxwiegand/predictosaurus/commit/d0919f2397564b90a66051154097957041157fc8))
* Calculate scores for haplotypes ([38be031](https://github.com/fxwiegand/predictosaurus/commit/38be0312f987f5844f87ff925aa371dbfc5fcafb))
* Calculate scores for haplotypes ([12518fa](https://github.com/fxwiegand/predictosaurus/commit/12518fa1bd7a35677436836a1219f47324caacba))
* Calculate stop penalty ([0aa1914](https://github.com/fxwiegand/predictosaurus/commit/0aa19147fb9300892c7403493aed2afe965420c1))
* Finish basic score implementation and restructure part of code ([2ba3d76](https://github.com/fxwiegand/predictosaurus/commit/2ba3d76259076d7428422e99a2924ce2f68f4636))
* Implement Epstein distance matrix and lookup for amino acids ([513c0db](https://github.com/fxwiegand/predictosaurus/commit/513c0db004e143b9f9ecc7571dee1236b661c1ae))
* Implement Miyata distance metric and add tests ([4002801](https://github.com/fxwiegand/predictosaurus/commit/40028018d6af78e02773ff5e8a12889b4d16c5cf))
* Start working on actual score calculation (WIP) ([24a1262](https://github.com/fxwiegand/predictosaurus/commit/24a1262d419169e8fbe283abc62d7e22f3625c02))
* Update phase according to haplotype properly considering coding ([3ae5746](https://github.com/fxwiegand/predictosaurus/commit/3ae5746e7d5a9e08748469c19bbdf3544a67237e))
* Use amino acid distance in score calculation ([70edff7](https://github.com/fxwiegand/predictosaurus/commit/70edff790032fdb75d5652703a343914ecfd4d2f))


### Bug Fixes

* Fix comment name ([0ca41b4](https://github.com/fxwiegand/predictosaurus/commit/0ca41b4445feb599fc0e662242ffbc06af2f4b92))
* Handle stop codons in amino acid distance metrics ([b0913bf](https://github.com/fxwiegand/predictosaurus/commit/b0913bf41dde2ac2c5983a91af7d18a46664e640))
* Penalize stop codons instead of methionine in scoring ([3cf8c29](https://github.com/fxwiegand/predictosaurus/commit/3cf8c29d197b44908231e95dba5683937ed20bc2))
* Remove unused import of std::ptr::write ([8307d35](https://github.com/fxwiegand/predictosaurus/commit/8307d35b227af363dbb030be44b440e5ef540144))

## [0.3.2](https://github.com/fxwiegand/predictosaurus/compare/v0.3.1...v0.3.2) (2025-04-29)


### Bug Fixes

* Remove nodes_in_between check for creating edges ([f1535ad](https://github.com/fxwiegand/predictosaurus/commit/f1535ad857f65319d63a47ad975c5dad6095d580))

## [0.3.1](https://github.com/fxwiegand/predictosaurus/compare/v0.3.0...v0.3.1) (2025-04-29)


### Bug Fixes

* Filter unique paths ([6106aba](https://github.com/fxwiegand/predictosaurus/commit/6106abaf93b85f60d05af801e868573f28f5af5b))

## [0.3.0](https://github.com/fxwiegand/predictosaurus/compare/v0.2.10...v0.3.0) (2025-04-23)


### Features

* Add subcommand for neoantigen prediction via peptides ([#78](https://github.com/fxwiegand/predictosaurus/issues/78)) ([4523e4f](https://github.com/fxwiegand/predictosaurus/commit/4523e4f0026fb616f3ea51bd37bc311efc02818a))


### Bug Fixes

* Cap numerical overshoot when calculating prob_present() ([95180ec](https://github.com/fxwiegand/predictosaurus/commit/95180ecabe51fa4e935d09cb2fec1c92f8d796b7))
* Cap numerical overshoot when calculating prob_present() ([43eefbe](https://github.com/fxwiegand/predictosaurus/commit/43eefbed4014892cad7e3bb69defde516499c815))
* Cap numerical overshoot when calculating prob_present() ([7f2e73c](https://github.com/fxwiegand/predictosaurus/commit/7f2e73c85c8653af2f19a592fd1cf1dc839e8f2f))
* Fix path generation ([862dd5d](https://github.com/fxwiegand/predictosaurus/commit/862dd5dd453185e561bcf7c964e334dadba58ac6))
* Fix path generation ([f440d26](https://github.com/fxwiegand/predictosaurus/commit/f440d262413ffea7022e2ab7c570842986674207))
* Fix weights concatenation ([31d9a12](https://github.com/fxwiegand/predictosaurus/commit/31d9a12fe03d7fd300295b229f502711027fde34))


### Performance Improvements

* Enhance transcript processing with variant detection ([#80](https://github.com/fxwiegand/predictosaurus/issues/80)) ([413db7f](https://github.com/fxwiegand/predictosaurus/commit/413db7f8697b161da019cdfdeae7f6b80525952f))

## [0.2.10](https://github.com/fxwiegand/predictosaurus/compare/v0.2.9...v0.2.10) (2025-02-14)


### Bug Fixes

* Consider splicing for each transcript by joining paths for CDSs ([#72](https://github.com/fxwiegand/predictosaurus/issues/72)) ([5c766c9](https://github.com/fxwiegand/predictosaurus/commit/5c766c9b1c8b0466da46961dad82d4efd5cc5ac2))

## [0.2.9](https://github.com/fxwiegand/predictosaurus/compare/v0.2.8...v0.2.9) (2025-01-17)


### Bug Fixes

* Adjust predictosaurus for non-unique CDS identifiers ([#60](https://github.com/fxwiegand/predictosaurus/issues/60)) ([e9c54ed](https://github.com/fxwiegand/predictosaurus/commit/e9c54ed58b4d33dd978b5d6ea50d8670da441bed))
* Allow  Ns in reverse_complement ([c684ece](https://github.com/fxwiegand/predictosaurus/commit/c684ece0c47cec7ff0a3f7435c73611bf1ee00e6))
* Remove edges between nodes with same position ([8879214](https://github.com/fxwiegand/predictosaurus/commit/8879214419cbb4299da7042e0593228e8b38853e))
* Remove edges between nodes with same position. ([1ab8ede](https://github.com/fxwiegand/predictosaurus/commit/1ab8edecb05bfb90f035ab5d6ac6392bc9b9acb8))
* Remove unused and potentially breaking placeholder nodes ([1d6f6bf](https://github.com/fxwiegand/predictosaurus/commit/1d6f6bf9377086b180cd4d6cf5f22789c7b6d468))
* Reverse reference sequence as expected by amino acid translations ([590c415](https://github.com/fxwiegand/predictosaurus/commit/590c4159498b6276f801af2e42159601faee2bde))

## [0.2.8](https://github.com/fxwiegand/predictosaurus/compare/v0.2.7...v0.2.8) (2025-01-07)


### Miscellaneous Chores

* release 0.2.8 ([1668725](https://github.com/fxwiegand/predictosaurus/commit/1668725766ee90157cb67c3d4c1a30980d0c66bb))

## [0.2.7](https://github.com/fxwiegand/predictosaurus/compare/v0.2.6...v0.2.7) (2025-01-07)


### Bug Fixes

* Ensure `vaf` scale domain is fixed between 0 and 1 ([01fd353](https://github.com/fxwiegand/predictosaurus/commit/01fd353a2a86064a9c5d939b1ea2c979e36a1d2a))

## [0.2.6](https://github.com/fxwiegand/predictosaurus/compare/v0.2.5...v0.2.6) (2024-12-16)


### Bug Fixes

* Create output dir if non existent ([79b1383](https://github.com/fxwiegand/predictosaurus/commit/79b1383075071d506bbfbeff13e2622a5295dd54))
* Fix feature graph generation ([283feec](https://github.com/fxwiegand/predictosaurus/commit/283feec842e031a288dacbdbd1e28fd0a4ebfb5d))

## [0.2.5](https://github.com/fxwiegand/predictosaurus/compare/v0.2.4...v0.2.5) (2024-12-12)


### Bug Fixes

* Remove primary key constraint on nodes table ([e62ac37](https://github.com/fxwiegand/predictosaurus/commit/e62ac3759bfbc9ae410b6101b82f53485d28be1b))

## [0.2.4](https://github.com/fxwiegand/predictosaurus/compare/v0.2.3...v0.2.4) (2024-12-12)


### Bug Fixes

* Fix wrong interpretation of probabilities ([7aeadb4](https://github.com/fxwiegand/predictosaurus/commit/7aeadb4a81462c6b1afd8ec8fba5506d39ac8f23))


### Performance Improvements

* Add indexes for write_graphs tables ([731ea91](https://github.com/fxwiegand/predictosaurus/commit/731ea91120e6ea8aee835aa0a01e02207c566cbf))

## [0.2.3](https://github.com/fxwiegand/predictosaurus/compare/v0.2.2...v0.2.3) (2024-12-10)


### Bug Fixes

* Fix multiple output file creation ([5ec7a40](https://github.com/fxwiegand/predictosaurus/commit/5ec7a40d6faef16745523d309927aa2f374c29e7))

## [0.2.2](https://github.com/fxwiegand/predictosaurus/compare/v0.2.1...v0.2.2) (2024-12-06)


### Bug Fixes

* Improve error messages ([c00a01f](https://github.com/fxwiegand/predictosaurus/commit/c00a01f64d205d53c8837b8381da758c8552199b))
* Improve error messages ([3f41035](https://github.com/fxwiegand/predictosaurus/commit/3f41035a54c39d093eea1c6ee52f2362ec8006df))

## [0.2.1](https://github.com/fxwiegand/predictosaurus/compare/v0.2.0...v0.2.1) (2024-12-06)


### Bug Fixes

* Fix out dir creation ([19dad8d](https://github.com/fxwiegand/predictosaurus/commit/19dad8dfec95fc7d742f645551a4b4c9860b60a3))


### Miscellaneous Chores

* release 0.2.1 ([054d16a](https://github.com/fxwiegand/predictosaurus/commit/054d16af036fb56ec93bb016ca99cc8cc84126fb))

## [0.2.0](https://github.com/fxwiegand/predictosaurus/compare/v0.1.1...v0.2.0) (2024-12-05)


### Features

* Add `--min-prop-present` parameter to build subcommand ([c4542cb](https://github.com/fxwiegand/predictosaurus/commit/c4542cb108d81db203ef1fbb0ab2f55b452ec97f))


### Bug Fixes

* Ignore samples present in calls file that are not present in given observations ([2bff91a](https://github.com/fxwiegand/predictosaurus/commit/2bff91a88314031c99d811b46e87d429da52702f))

## [0.1.1](https://github.com/fxwiegand/predictosaurus/compare/v0.1.0...v0.1.1) (2024-12-04)


### Miscellaneous Chores

* release 0.1.1 ([8551a11](https://github.com/fxwiegand/predictosaurus/commit/8551a110d543be8b84acdfe0f0e0e2692f17ec15))

## 0.1.0 (2024-12-04)


### Features

* Add reverse paths for handling features on reverse strand ([897bd81](https://github.com/fxwiegand/predictosaurus/commit/897bd81193fd062b7adae954e373660ecb511607))
* Add support for features on backward strand ([e6b2141](https://github.com/fxwiegand/predictosaurus/commit/e6b21411d3fdce1b54a36f115e8cb76d099ea030))
* Allow conversion of haplotype path to protein ([3fa2c75](https://github.com/fxwiegand/predictosaurus/commit/3fa2c7598268cf7f197aa8943f25b992fb80780e))
* Calculate 0-based position of variant on reverse strand ([8321619](https://github.com/fxwiegand/predictosaurus/commit/8321619843be637c57aad43d671701a93d85b971))
* Enable serialization and deserialization for Graph ([200b4dd](https://github.com/fxwiegand/predictosaurus/commit/200b4dd8cba77efcd3aef66c0151d0bc845c7488))
* Implement show command ([ba4ea1d](https://github.com/fxwiegand/predictosaurus/commit/ba4ea1d1776fdab0d180a9cb96d5b1c62fca9122))
* Implement show command ([a2f243e](https://github.com/fxwiegand/predictosaurus/commit/a2f243ecca5cafd9f65724e103acece84b4bfac8))
* Implement show command ([3bc8d30](https://github.com/fxwiegand/predictosaurus/commit/3bc8d3080dae4e5106a9b7113e08fa92c514def3))
* Implement show command for HTML format ([4af776b](https://github.com/fxwiegand/predictosaurus/commit/4af776b69a00b41a73ec8d8e874195d61be291cd))
* Implement show command for tsv format ([87be0d9](https://github.com/fxwiegand/predictosaurus/commit/87be0d984711520da761e534b540fca3f7a3b8f6))
* Implement show command for tsv format ([4616bc6](https://github.com/fxwiegand/predictosaurus/commit/4616bc6b70c4df1916f8ff68301fe5f2756ee249))
* Implement show command for tsv format ([f0a8965](https://github.com/fxwiegand/predictosaurus/commit/f0a89655851eaf15628d247d6a5857c825e4f87a))
* Restructure into multiple subcommands ([#13](https://github.com/fxwiegand/predictosaurus/issues/13)) ([4462a23](https://github.com/fxwiegand/predictosaurus/commit/4462a233f037cdccc87ec5acfc40af5f0f69d458))
* Support insertions longer 3 ([252d1d6](https://github.com/fxwiegand/predictosaurus/commit/252d1d64577ab4046b0edc4a2795ba3e2d970932))
* Test graph serialization using serde_json ([5074b49](https://github.com/fxwiegand/predictosaurus/commit/5074b49001a8228c67a3c645ae73469812823cec))
* Use duckdb to write graphs ([#35](https://github.com/fxwiegand/predictosaurus/issues/35)) ([19955ea](https://github.com/fxwiegand/predictosaurus/commit/19955eadaac9d3c0f9dadb7480e8265dca71a285))


### Bug Fixes

* Create path table in bd before writing to it ([c335418](https://github.com/fxwiegand/predictosaurus/commit/c3354189f046ce9a24cd78baa5b5349435a790fe))
* Fix impact color scheme ([b8bed5e](https://github.com/fxwiegand/predictosaurus/commit/b8bed5e733d49863efead1100ad1bb18c15755b9))
* Fix impact for amino acid altering variants ([38548a0](https://github.com/fxwiegand/predictosaurus/commit/38548a0e35ecfa95db48069a0c36012b61365abb))
* Handle variants at start and end of reference ([794a7e0](https://github.com/fxwiegand/predictosaurus/commit/794a7e081eb1775cb1b9c83639c080980efda6b3))
* Refactor field names ([2378127](https://github.com/fxwiegand/predictosaurus/commit/2378127b45be953b96f2f61aef799bf6a8e2f6fa))
* Refine behaviour for invalid EventProbs ([c1a300d](https://github.com/fxwiegand/predictosaurus/commit/c1a300d7ba7f4462e2774a869e31befd0d945c3b))
* Refine behaviour for invalid EventProbs ([4342223](https://github.com/fxwiegand/predictosaurus/commit/4342223da59defd6eb39c76fce4225765dc3efa9))
* Refine behaviour for invalid EventProbs ([ac541b1](https://github.com/fxwiegand/predictosaurus/commit/ac541b1382f5304913c8aae762e663a96e63151d))
* Refine behaviour for invalid EventProbs ([d06b7c4](https://github.com/fxwiegand/predictosaurus/commit/d06b7c45849c604f7dbb9704ddc89290b9c08620))
* Refine behaviour for invalid EventProbs ([89ae172](https://github.com/fxwiegand/predictosaurus/commit/89ae172f7948ff94469d9a0f7dac4b6615690f72))
* Remove duplicate ref nodes for multiple variants with same position ([daf9ee0](https://github.com/fxwiegand/predictosaurus/commit/daf9ee0718050cd4cf9782471f82fa036fd05408))
* Use json5 to serialize and deserialize EventProbs ([3330305](https://github.com/fxwiegand/predictosaurus/commit/333030533326cb24acaec31ff776014a60b7cbc9))
