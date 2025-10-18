(function () {
  "use strict";

  function ready(callback) {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", callback);
    } else {
      callback();
    }
  }

  function toArray(nodeList) {
    return Array.prototype.slice.call(nodeList || []);
  }

  function parseTargets(section) {
    const raw = section.getAttribute("data-platform-section") || "";
    if (!raw.trim()) {
      return [];
    }

    return raw
      .split(/[,\s]+/)
      .map(function (item) {
        return item.trim();
      })
      .filter(Boolean);
  }

  function toggleSection(section, isActive) {
    section.classList.toggle("d-none", !isActive);

    const controls = toArray(
      section.querySelectorAll("input, select, textarea, button")
    );
    controls.forEach(function (control) {
      control.disabled = !isActive;
    });
  }

  function applySelection(select, sections) {
    const selected = select.value;

    sections.forEach(function (section) {
      const targets = parseTargets(section);
      const matches =
        targets.length === 0 || targets.indexOf(selected) !== -1;

      toggleSection(section, matches);
    });
  }

  ready(function () {
    const select = document.querySelector("#id_sequencing_platform");
    if (!select) {
      return;
    }

    const sections = toArray(
      document.querySelectorAll("[data-platform-section]")
    );
    if (sections.length === 0) {
      return;
    }

    const handler = function () {
      applySelection(select, sections);
    };

    select.addEventListener("change", handler);
    select.addEventListener("input", handler);

    handler();
  });
})();
