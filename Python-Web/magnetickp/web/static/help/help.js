const searchInput = document.getElementById("docSearch");
const noResults = document.getElementById("noResults");
const navToggle = document.getElementById("navToggle");
const sidebar = document.getElementById("sidebar");
const languageToggle = document.getElementById("languageToggle");
const sections = Array.from(document.querySelectorAll(".doc-section"));
const navLinks = Array.from(document.querySelectorAll(".sidebar a[href^='#']"));

let currentLanguage = "zh";
try {
  currentLanguage = localStorage.getItem("magnetickp-language") === "en" ? "en" : "zh";
} catch (_) {}

function renderCodeComments() {
  document.querySelectorAll("code.language-python").forEach((code) => {
    const source = code.dataset.bilingualSource || code.textContent;
    if (!code.dataset.bilingualSource) code.dataset.bilingualSource = source;
    code.textContent = source.split("\n").map((line) => {
      const separatorIndex = line.indexOf(" || ");
      if (separatorIndex === -1) return line;
      if (currentLanguage === "zh") return line.slice(0, separatorIndex);
      const indentation = line.match(/^\s*/)?.[0] || "";
      return `${indentation}# ${line.slice(separatorIndex + 4)}`;
    }).join("\n");
  });
}

function setLanguage(language) {
  currentLanguage = language === "en" ? "en" : "zh";
  document.body.dataset.language = currentLanguage;
  document.documentElement.lang = currentLanguage === "zh" ? "zh-CN" : "en";
  document.title = currentLanguage === "zh"
    ? "MagneticKP · Python + Web 帮助"
    : "MagneticKP · Python + Web Documentation";
  searchInput.placeholder = currentLanguage === "zh"
    ? "搜索教程、函数或示例…"
    : "Search tutorials, functions, or examples…";
  languageToggle.textContent = currentLanguage === "zh" ? "EN" : "中文";
  languageToggle.setAttribute(
    "aria-label",
    currentLanguage === "zh" ? "Switch to English" : "切换至中文",
  );
  renderCodeComments();
  try { localStorage.setItem("magnetickp-language", currentLanguage); } catch (_) {}
}

function addCopyButtons() {
  document.querySelectorAll("pre").forEach((block) => {
    if (block.querySelector(".copy-button")) return;
    const button = document.createElement("button");
    button.type = "button";
    button.className = "copy-button";
    button.textContent = currentLanguage === "zh" ? "复制" : "Copy";
    button.addEventListener("click", async () => {
      const code = block.querySelector("code");
      if (!code) return;
      await navigator.clipboard.writeText(code.textContent);
      button.textContent = currentLanguage === "zh" ? "已复制" : "Copied";
      button.classList.add("copied");
      window.setTimeout(() => {
        button.textContent = currentLanguage === "zh" ? "复制" : "Copy";
        button.classList.remove("copied");
      }, 1200);
    });
    block.append(button);
  });
}

function filterSections() {
  const query = searchInput.value.trim().toLocaleLowerCase();
  let visibleCount = 0;
  sections.forEach((section) => {
    const visible = !query || section.textContent.toLocaleLowerCase().includes(query);
    section.hidden = !visible;
    if (visible) visibleCount += 1;
    const link = navLinks.find((item) => item.hash === `#${section.id}`);
    if (link) link.hidden = !visible;
  });
  noResults.hidden = visibleCount !== 0;
}

function setActiveLink(id) {
  navLinks.forEach((link) => link.classList.toggle("active", link.hash === `#${id}`));
}

const observer = "IntersectionObserver" in window
  ? new IntersectionObserver((entries) => {
      const visible = entries
        .filter((entry) => entry.isIntersecting && !entry.target.hidden)
        .sort((left, right) => left.boundingClientRect.top - right.boundingClientRect.top);
      if (visible[0]) setActiveLink(visible[0].target.id);
    }, { rootMargin: "-18% 0px -68% 0px", threshold: 0 })
  : null;

sections.forEach((section) => observer?.observe(section));
searchInput.addEventListener("input", filterSections);
searchInput.addEventListener("keydown", (event) => {
  if (event.key !== "Escape") return;
  searchInput.value = "";
  filterSections();
  searchInput.blur();
});

document.addEventListener("keydown", (event) => {
  if (event.key === "/" && document.activeElement !== searchInput) {
    event.preventDefault();
    searchInput.focus();
  }
});

navToggle.addEventListener("click", () => {
  const expanded = navToggle.getAttribute("aria-expanded") === "true";
  navToggle.setAttribute("aria-expanded", String(!expanded));
  sidebar.classList.toggle("open", !expanded);
});

navLinks.forEach((link) => link.addEventListener("click", () => {
  sidebar.classList.remove("open");
  navToggle.setAttribute("aria-expanded", "false");
}));

languageToggle.addEventListener("click", () => {
  setLanguage(currentLanguage === "zh" ? "en" : "zh");
  document.querySelectorAll(".copy-button").forEach((button) => {
    if (!button.classList.contains("copied")) {
      button.textContent = currentLanguage === "zh" ? "复制" : "Copy";
    }
  });
});

setLanguage(currentLanguage);
addCopyButtons();
